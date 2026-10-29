"""Run a command repeatedly and measure it: wall and CPU time, peak memory.

Each run is a child process waited on with os.wait4, so its CPU time is its
own (and its children's, for a shell or an R script that runs a binary); its
peak resident set size comes from GNU time.
"""

from __future__ import annotations

import contextlib
import os
import platform
import statistics
import subprocess
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path

# GNU time (`time -f %M` gives the peak RSS in kilobytes); Linux
GNU_TIME = os.environ.get("POPMIX_GNU_TIME", "/usr/bin/time")


@dataclass
class Run:
    wall_s: float
    cpu_s: float
    peak_rss_mb: float


@dataclass
class Measured:
    runs: list[Run] = field(default_factory=list)

    def summary(self) -> dict:
        out = {"runs": [r.__dict__ for r in self.runs]}
        for key in ("wall_s", "cpu_s", "peak_rss_mb"):
            values = sorted(getattr(r, key) for r in self.runs)
            median = statistics.median(values)
            q = statistics.quantiles(values, n=4) if len(values) >= 2 else [median, median, median]
            mean = statistics.fmean(values)
            sd = statistics.stdev(values) if len(values) >= 2 else 0.0
            out[key] = {"median": median, "iqr": q[2] - q[0], "min": values[0], "max": values[-1],
                        "cv": sd / mean if mean else 0.0}
        return out


def run_once(cmd: list[str], cwd: Path, env: dict[str, str] | None = None, log: Path | None = None) -> Run:
    """Run cmd to completion; raise if it fails."""
    # GNU time reports the command's own peak RSS: a child forked from this
    # (large) Python process would start with this process's resident set
    with tempfile.NamedTemporaryFile("r", encoding="utf-8", suffix=".rss") as rss, contextlib.ExitStack() as stack:
        out = stack.enter_context(open(log, "ab")) if log else subprocess.DEVNULL
        start = time.perf_counter()
        proc = subprocess.Popen([GNU_TIME, "-f", "%M", "-o", rss.name, *cmd], cwd=cwd, env=env, stdout=out,
                                stderr=subprocess.STDOUT)
        _, status, usage = os.wait4(proc.pid, 0)
        wall = time.perf_counter() - start
        proc.returncode = os.waitstatus_to_exitcode(status)
        if proc.returncode != 0:
            raise RuntimeError(f"{' '.join(cmd)} exited {proc.returncode} (log: {log})")
        peak_kb = int(rss.read().split()[-1])
    return Run(wall_s=wall, cpu_s=usage.ru_utime + usage.ru_stime, peak_rss_mb=peak_kb / 1024)


def measure(cmd: list[str], cwd: Path, repeats: int, env: dict[str, str] | None = None,
            log: Path | None = None, warmup: bool = True) -> Measured:
    """Run cmd `repeats` times (after one unmeasured run that warms the page
    cache, so every repeat reads from memory)."""
    if warmup:
        run_once(cmd, cwd, env, log)
    m = Measured()
    for _ in range(repeats):
        m.runs.append(run_once(cmd, cwd, env, log))
    return m


def machine() -> dict:
    info = {"python": platform.python_version(), "system": platform.platform(), "cpus": os.cpu_count()}
    try:
        for line in Path("/proc/cpuinfo").read_text(encoding="utf-8").splitlines():
            if line.startswith("model name"):
                info["cpu"] = line.split(":", 1)[1].strip()
                break
        mem = Path("/proc/meminfo").read_text(encoding="utf-8").split()[1]
        info["memory_gb"] = round(int(mem) / 1024 / 1024, 1)
    except OSError:
        pass
    return info
