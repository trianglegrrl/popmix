"""The process rules in CLAUDE.md that a machine can check (TASKS.md K7, I7).

- Tests use real result objects, never mock result classes.
- An assertion states one expectation: no ``assert a or b``, which passes
  when either half does and so tests neither.
- Source files stay small: none past ``MAX_LINES``. The files over it when
  the rule came in are listed with their size then; they may only shrink,
  and leave the list once under the limit.
"""

from __future__ import annotations

import ast
from pathlib import Path

ROOT = Path(__file__).parents[2]
TEST_DIRS = [ROOT / "tests"]
SOURCE_DIRS = [ROOT / d for d in ("crates", "python", "tests", "benchmarks")]
MAX_LINES = 800

# files over MAX_LINES allowed to stay so, with their sizes, which may only
# go down: none since TASKS.md I7 split the 13 there were (2026-10-27)
TOO_LONG: dict[str, int] = {}


def _files(dirs: list[Path], suffixes: tuple[str, ...]) -> list[Path]:
    return sorted(p for d in dirs for p in d.rglob("*")
                  if p.suffix in suffixes and p.is_file()
                  and not {"target", "__pycache__", ".hypothesis"} & set(p.parts))


def _test_modules() -> list[tuple[Path, ast.Module]]:
    return [(p, ast.parse(p.read_text(), str(p))) for p in _files(TEST_DIRS, (".py",))]


def test_no_assertion_is_an_or_of_expectations():
    found = [f"{p.relative_to(ROOT)}:{node.lineno}"
             for p, tree in _test_modules() for node in ast.walk(tree)
             if isinstance(node, ast.Assert) and isinstance(node.test, ast.BoolOp)
             and isinstance(node.test.op, ast.Or)]
    assert found == []


def test_no_mock_result_classes():
    found = []
    for p, tree in _test_modules():
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef) and node.name.startswith(("Mock", "Fake")):
                found.append(f"{p.relative_to(ROOT)}:{node.lineno} class {node.name}")
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                names = [node.module or ""] if isinstance(node, ast.ImportFrom) else [a.name for a in node.names]
                if any(n == "unittest.mock" or n.startswith("unittest.mock.") or n == "mock" for n in names):
                    found.append(f"{p.relative_to(ROOT)}:{node.lineno} imports mock")
    assert found == []


def test_source_files_stay_small():
    sizes = {str(p.relative_to(ROOT)): len(p.read_text().splitlines())
             for p in _files(SOURCE_DIRS, (".rs", ".py"))}
    new = {f: n for f, n in sizes.items() if n > MAX_LINES and f not in TOO_LONG}
    grown = {f: (TOO_LONG[f], sizes[f]) for f in TOO_LONG if f in sizes and sizes[f] > TOO_LONG[f]}
    fixed = [f for f in TOO_LONG if sizes.get(f, 0) <= MAX_LINES]
    assert new == {}, f"split these, or keep them under {MAX_LINES} lines"
    assert grown == {}, "files over the limit may only shrink (was, now)"
    assert fixed == [], "under the limit (or gone) now: take them out of TOO_LONG"
