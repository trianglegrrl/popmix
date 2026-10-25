#!/usr/bin/env python3
"""Reference for popmix's parsers of classic ADMIXTOOLS output.

Runs ADMIXTOOLS 7.0.2 (build_classic.sh) qpDstat, qp3Pop, qpAdm, qpWave,
qpGraph, qpF4ratio and qpfstats on dg_6pop populations of AADR v62, keeps
their output in tests/golden/admixtools_classic/v7.0.2/, and records what
admixtools2's parsers (parse_qpdstat_output, parse_qp3pop_output,
parse_qpadm_output, parse_qpgraph_output, parse_qpgraph_graphfile,
parse_qpf4ratio_output, parse_qpfstats_output) make of each, in parsed.json.

Usage: python tests/golden/admixtools_classic/generate_classic.py
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

HERE = Path(__file__).parent
OUT = HERE / "v7.0.2"
BIN = Path(os.environ.get("ADMIXTOOLS_BIN", Path.home() / ".local/share/popmix-references/admixtools-7.0.2/bin"))
AADR = "/references/aadr/v62.0_1240k_public"
POPS = ["Chimp.REF", "Mbuti.DG", "Yoruba.DG", "French.DG", "Han.DG", "Papuan.DG"]
# the classic programs read the whole genotype file: run them on the
# populations they need, cut out with convertf (not kept in the repository)
SUBSET = Path(os.environ.get("POPMIX_REFERENCE_CACHE", Path.home() / ".cache/popmix-references")) / "aadr_dg6pop"

QPDSTAT = ["Chimp.REF Mbuti.DG French.DG Han.DG", "Chimp.REF Mbuti.DG Yoruba.DG French.DG",
           "Mbuti.DG Yoruba.DG Papuan.DG Han.DG"]
QP3POP = ["Han.DG Yoruba.DG French.DG", "Mbuti.DG Han.DG Papuan.DG"]
LEFT = ["French.DG", "Han.DG", "Yoruba.DG"]
RIGHT = ["Mbuti.DG", "Papuan.DG", "Chimp.REF"]
F4RATIO = ["Chimp.REF Mbuti.DG : French.DG Han.DG :: Chimp.REF Mbuti.DG : Yoruba.DG Han.DG"]
GRAPH = """root R
label Chimp.REF Chimp
label Mbuti.DG Mbuti
label Yoruba.DG Yoruba
label French.DG French
label Han.DG Han
edge a R Chimp
edge b R n1
edge c n1 Mbuti
edge d n1 n2
edge e n2 Yoruba
edge f n2 n3
edge g n3 eW
edge h n3 Han
admix adm eW n2
edge i adm French
"""


def par(**kw) -> str:
    base = {"genotypename": f"{SUBSET}.geno", "snpname": f"{SUBSET}.snp", "indivname": f"{SUBSET}.ind"}
    return "".join(f"{k}: {v}\n" for k, v in {**base, **kw}.items())


def run(prog: str, name: str, parfile: str, *args: str) -> Path:
    p = OUT / f"{name}.par"
    p.write_text(parfile)
    out = OUT / f"{name}.out"
    with open(out, "w") as fh:
        subprocess.run([str(BIN / prog), "-p", str(p), *args], stdout=fh, stderr=subprocess.STDOUT, check=True,
                       cwd=OUT)
    print("ran", prog, "->", out.name)
    return out


def subset():
    if Path(f"{SUBSET}.geno").exists():
        return
    SUBSET.parent.mkdir(parents=True, exist_ok=True)
    pops = SUBSET.parent / "aadr_dg6pop.pops"
    pops.write_text("\n".join(POPS) + "\n")
    p = SUBSET.parent / "aadr_dg6pop.par"
    p.write_text(f"genotypename: {AADR}.geno\nsnpname: {AADR}.snp\nindivname: {AADR}.ind\n"
                 f"outputformat: PACKEDANCESTRYMAP\ngenooutfilename: {SUBSET}.geno\n"
                 f"snpoutfilename: {SUBSET}.snp\nindoutfilename: {SUBSET}.ind\npoplistname: {pops}\n")
    subprocess.run([str(BIN / "convertf"), "-p", str(p)], check=True, stdout=subprocess.DEVNULL)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    subset()
    (OUT / "qpdstat.pops").write_text("\n".join(QPDSTAT) + "\n")
    (OUT / "qp3pop.pops").write_text("\n".join(QP3POP) + "\n")
    (OUT / "left.pops").write_text("\n".join(LEFT) + "\n")
    (OUT / "right.pops").write_text("\n".join(RIGHT) + "\n")
    (OUT / "f4ratio.pops").write_text("\n".join(F4RATIO) + "\n")
    (OUT / "graph.txt").write_text(GRAPH)
    # qpfstats needs heterozygous sites in its base (first) population; Chimp.REF has none
    (OUT / "fstats.pops").write_text("\n".join(["Mbuti.DG", "Yoruba.DG", "French.DG", "Han.DG", "Chimp.REF"]) + "\n")

    run("qpDstat", "qpdstat", par(popfilename="qpdstat.pops", f4mode="YES"))
    run("qp3Pop", "qp3pop", par(popfilename="qp3pop.pops"))
    run("qpAdm", "qpadm", par(popleft="left.pops", popright="right.pops", details="YES"))
    run("qpWave", "qpwave", par(popleft="left.pops", popright="right.pops"))
    run("qpF4ratio", "qpf4ratio", par(popfilename="f4ratio.pops"))
    run("qpGraph", "qpgraph", par(outpop="NULL", useallsnps="NO", blgsize="0.05", forcezmode="YES",
                                  lsqmode="NO", diag=".0001", bigiter="6", hires="YES", lambdascale="1"),
        "-g", "graph.txt", "-o", "qpgraph.ggg", "-d", "qpgraph.dot")
    run("qpfstats", "qpfstats", par(poplistname="fstats.pops", fstatsoutname="qpfstats.fstats"))

    r = f"""
    suppressPackageStartupMessages({{library(admixtools); library(jsonlite)}})
    setwd("{OUT}")
    tab <- function(d) lapply(split(as.data.frame(d), seq_len(nrow(d))), as.list)
    num <- function(d) {{ d <- as.data.frame(d); for (k in names(d)) if (is.numeric(d[[k]])) d[[k]] <- sprintf("%.17g", d[[k]]); unname(tab(d)) }}
    qa <- admixtools:::parse_qpadm_output("qpadm.out")
    qg <- admixtools:::parse_qpgraph_output("qpgraph.ggg", logfile = "qpgraph.out")
    qf <- admixtools:::parse_qpfstats_output("qpfstats.out")
    out <- list(
      qpdstat = num(admixtools:::parse_qpdstat_output("qpdstat.out")),
      qp3pop = num(admixtools:::parse_qp3pop_output("qp3pop.out")),
      qpadm = list(weights = num(qa$weights), popdrop = num(qa$popdrop)),
      qpf4ratio = num(admixtools:::parse_qpf4ratio_output("qpf4ratio.out")),
      qpgraph = list(edges = num(qg$edges), score = sprintf("%.17g", qg$score), f2 = num(qg$f2), f3 = num(qg$f3),
                     f4 = if (is.null(qg$f4)) NULL else num(qg$f4)),
      qpgraph_edges_ggg = num(admixtools:::parse_qpgraph_output_edges("qpgraph.ggg")),
      # split_multi = TRUE (the default) fails in 2.0.10 (DEVIATIONS.md, ADMIXTOOLS2-0005)
      graphfile = num(admixtools:::parse_qpgraph_graphfile("graph.txt", split_multi = FALSE)),
      qpfstats = list(fstat = if (is.null(qf$fstat)) NULL else num(qf$fstat), basis = num(qf$basis),
                      jest = if (is.null(qf$jest)) NULL else num(qf$jest)))
    writeLines(toJSON(out, auto_unbox = TRUE, null = "null", na = "null"), "parsed.json")
    """
    subprocess.run(["/usr/bin/Rscript", "-e", r], check=True)
    print("wrote", OUT / "parsed.json")


if __name__ == "__main__":
    main()
