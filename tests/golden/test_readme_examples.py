"""The README's example outputs are what the commands print (TASKS.md M5).

Each ```bash block followed by an **Output** block holds one command; run
in a directory where ``aadr_v66p1`` is AADR v66.p1 labelled by data type
(tests/golden/aadr/build_aadr_v66p1.py) and ``my_samples``/``reference``
are tests/fixtures/merge_a and merge_b, it prints that output, exactly. A
``...`` line in the output stands for any run of lines. An output marked
"on stderr" is the command's error, with a nonzero exit.
"""

from __future__ import annotations

import os
import re
import shlex
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[2]
AADR = Path(os.environ.get("POPMIX_REFERENCE_CACHE", Path.home() / ".cache/popmix-references")) / "aadr_v66p1" / "aadr_v66p1"
FIXTURES = ROOT / "tests" / "fixtures"
PAIR = re.compile(r"```bash\n([^`]*)```\n\n(\*\*Output[^\n]*)\n```[a-z]*\n([^`]*)```")


def _examples() -> list[tuple[str, str, str]]:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    return [(code, header, output) for code, header, output in PAIR.findall(readme)]


def _command(code: str) -> str:
    lines = [line for line in code.replace("\\\n", " ").splitlines() if line.strip() and not line.startswith("#")]
    assert len(lines) == 1, f"an example with an output runs one command: {lines}"
    return lines[0]


def _matches(expected: list[str], actual: list[str]) -> bool:
    """Whether ``actual`` is ``expected``, each "..." line standing for any lines."""
    parts: list[list[str]] = [[]]
    for line in expected:
        if line == "...":
            parts.append([])
        else:
            parts[-1].append(line)
    if len(parts) == 1:
        return actual == expected
    first, *middle, last = parts
    if actual[:len(first)] != first or len(actual) < len(first) + len(last):
        return False
    if last and actual[len(actual) - len(last):] != last:
        return False
    i, end = len(first), len(actual) - len(last)
    for part in middle:
        while i + len(part) <= end and actual[i:i + len(part)] != part:
            i += 1
        if i + len(part) > end:
            return False
        i += len(part)
    return True


def _lines(text: str) -> list[str]:
    return [line.rstrip() for line in text.rstrip("\n").split("\n")]


def test_the_readme_has_its_examples():
    commands = [_command(code).split()[1] for code, _, _ in _examples()]
    for name in ("info", "populations", "merge", "f2", "f3", "f4", "qpadm", "qpwave", "dates"):
        assert name in commands


def test_the_ellipsis_matches_any_run_of_lines():
    assert _matches(["a", "...", "d"], ["a", "b", "c", "d"])
    assert _matches(["a", "...", "d"], ["a", "d"])
    assert not _matches(["a", "...", "d"], ["a", "b"])
    assert _matches(["a", "...", "c", "...", "e"], ["a", "b", "c", "d", "e"])
    assert not _matches(["a", "b"], ["a", "b", "c"])


@pytest.fixture(scope="module")
def example_dir(tmp_path_factory) -> Path:
    if not Path(f"{AADR}.geno").exists():
        pytest.skip("AADR v66.p1 dataset not built")
    d = tmp_path_factory.mktemp("readme")
    for ext in ("geno", "snp", "ind"):
        (d / f"aadr_v66p1.{ext}").symlink_to(f"{AADR}.{ext}")
        (d / f"my_samples.{ext}").symlink_to(FIXTURES / f"merge_a.{ext}")
        (d / f"reference.{ext}").symlink_to(FIXTURES / f"merge_b.{ext}")
    return d


@pytest.mark.parametrize("code,header,output", _examples(), ids=lambda x: x.split("\n")[0][:60])
def test_the_example_prints_its_output(example_dir, code, header, output):
    run = subprocess.run(shlex.split(_command(code)), cwd=example_dir, capture_output=True, text=True)
    if "on stderr" in header:
        assert run.returncode != 0
        printed = run.stderr
    else:
        assert run.returncode == 0, run.stderr
        printed = run.stdout
    assert _matches(_lines(output), _lines(printed)), printed
