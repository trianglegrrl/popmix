"""Every `popmix ...` command line in README.md, docs/guide.md and TESTING.md
exists (D4).

Each documented line names a real command, and every flag it uses is one of
that command's options (or a global one). The documentation used to show
commands (qpadm-batch, dates-batch) and flags (--report) that do not exist.
"""

from __future__ import annotations

import re
import shlex
from pathlib import Path

import click
import pytest

from popmix.cli import cli

ROOT = Path(__file__).resolve().parents[2]
DOCS = ["README.md", "docs/guide.md", "TESTING.md"]


def _command_lines(path: Path) -> list[tuple[int, str]]:
    """`popmix ...` lines inside fenced shell blocks, continuations joined."""
    lines, out, in_block, lang = path.read_text(encoding="utf-8").splitlines(), [], False, ""
    i = 0
    while i < len(lines):
        line = lines[i]
        fence = re.match(r"^\s*```(\w*)", line)
        if fence:
            in_block, lang = (not in_block), fence.group(1)
            i += 1
            continue
        if in_block and lang in ("", "bash", "sh", "shell", "console"):
            start = i
            text = line.strip()
            while text.endswith("\\") and i + 1 < len(lines):
                i += 1
                text = text[:-1] + " " + lines[i].strip()
            text = re.sub(r"^\$\s+", "", text)
            if text.startswith("popmix "):
                out.append((start + 1, text.split(" #")[0].split(" | ")[0].strip()))
        i += 1
    return out


def _flags(command: click.Command) -> set[str]:
    return {opt for p in command.params for opt in (*p.opts, *p.secondary_opts)}


CASES = [(doc, n, text) for doc in DOCS if (ROOT / doc).exists() for n, text in _command_lines(ROOT / doc)]


@pytest.mark.parametrize(("doc", "line", "text"), CASES, ids=[f"{d}:{n}" for d, n, _ in CASES])
def test_documented_command_exists(doc, line, text):
    words = shlex.split(text.replace("...", "X").replace("…", "X"))[1:]
    global_flags = _flags(cli) | {"--help"}
    # global options before the command (--threads 8)
    while words and words[0].startswith("-"):
        flag = words.pop(0).split("=")[0]
        assert flag in global_flags, f"{doc}:{line}: unknown global option {flag}"
        if flag in ("--threads", "-t") and words:
            words.pop(0)
    if not words or words[0] in ("--help", "--version"):
        return
    name = words.pop(0)
    command = cli.commands.get(name)
    assert command is not None, f"{doc}:{line}: no command {name!r}"
    known = _flags(command) | {"--help"}
    for word in words:
        if word.startswith("-") and not re.fullmatch(r"-\d.*", word):
            flag = word.split("=")[0]
            assert flag in known, f"{doc}:{line}: {name} has no option {flag}"


def test_some_commands_are_documented():
    assert len(CASES) > 20
