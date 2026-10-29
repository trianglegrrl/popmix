"""The agent skill (TASKS.md N1-N4) documents the CLI exactly, and its examples run.

The skill is ``python/popmix/skill/`` (``skills/popmix`` in the repository):
``SKILL.md``, ``python-api.md`` (tests/python/test_skill_api.py) and one file per command,
``commands/<command>.md`` (a subcommand of a group would be
``<group>-<command>.md``). These tests walk the Click tree, so a command or
option added, renamed or removed without its documentation fails here:

- every command has exactly one file and every file one command;
- every argument and option has a row in its file's options table, under
  each of its flag spellings, with its type and its default;
- every ``popmix`` line in a shell block names a real command and real flags;
- every shell block is marked ``# runs`` or ``# not run: <reason>``, and the
  ``# runs`` blocks execute (bash, ``popmix`` on PATH) in a directory holding
  the small bundled data that SKILL.md describes, and exit 0;
- the skill states popmix's version, ships inside the package and
  ``popmix install-skill`` copies it.
"""

from __future__ import annotations

import os
import re
import shlex
import shutil
import subprocess
import sys
from importlib import resources
from pathlib import Path

import click
import pytest

import popmix as pm
from popmix.cli import cli

try:
    from click._utils import UNSET
except ImportError:  # click before 8.3: no default is None
    UNSET = None

ROOT = Path(__file__).resolve().parents[2]
SKILL = ROOT / "python" / "popmix" / "skill"
COMMANDS = SKILL / "commands"
FAMILY = ROOT / "tests" / "golden" / "dates" / "family_chr20_22"
FIXTURES = ROOT / "tests" / "fixtures"


def walk(group: click.Group, path: tuple[str, ...] = ()) -> dict[str, tuple[tuple[str, ...], click.Command]]:
    """Every leaf command, keyed by its file stem (the path joined by '-')."""
    out = {}
    for name, command in group.commands.items():
        if isinstance(command, click.Group):
            out.update(walk(command, (*path, name)))
        else:
            out["-".join((*path, name))] = ((*path, name), command)
    return out


ALL = {stem: command for stem, (_, command) in walk(cli).items()}
PATHS = {stem: path for stem, (path, _) in walk(cli).items()}


def doc(stem: str) -> str:
    return (COMMANDS / f"{stem}.md").read_text(encoding="utf-8")


# ---------------------------------------------------------------- tables


def table_rows(text: str) -> list[list[str]]:
    """Cells of every markdown table row (header and separator rows dropped)."""
    rows = []
    for line in text.splitlines():
        if not line.startswith("|"):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if all(re.fullmatch(r":?-+:?", c) for c in cells) or cells[0] in ("Option", "Argument"):
            continue
        rows.append(cells)
    return rows


def row_for(param: click.Parameter, rows: list[list[str]]) -> list[str] | None:
    names = spellings(param)
    for cells in rows:
        ticked = set(re.findall(r"`([^`]+)`", cells[0]))
        if names <= ticked:
            return cells
    return None


def spellings(param: click.Parameter) -> set[str]:
    if isinstance(param, click.Argument):
        return {param.human_readable_name}
    return {*param.opts, *param.secondary_opts}


def type_words(param: click.Parameter) -> list[str]:
    """What the type cell (or, for a choice, the row) must say."""
    t = param.type
    if isinstance(param, click.Option) and param.is_flag:
        return ["flag"]
    if isinstance(t, click.Choice):
        return list(t.choices)
    if isinstance(t, (click.types.IntParamType, click.IntRange)):
        return ["integer"]
    if isinstance(t, (click.types.FloatParamType, click.FloatRange)):
        return ["number"]
    if isinstance(t, click.Path):
        return ["path"]
    return []  # text: "text" or "path", both accepted below


def check_row(stem: str, param: click.Parameter, cells: list[str]) -> list[str]:
    errors = []
    label = f"{stem}: {'/'.join(sorted(spellings(param)))}"
    if len(cells) < 4:
        return [f"{label}: row needs Option | Type | Default | Meaning"]
    kind, default_cell, row = cells[1].lower(), cells[2], " ".join(cells)
    words = type_words(param)
    if isinstance(param.type, click.Choice):
        missing = [c for c in words if f"`{c}`" not in row]
        if missing:
            errors.append(f"{label}: choices {missing} not in its row")
    elif words and words[0] not in kind:
        errors.append(f"{label}: type cell {cells[1]!r} should say {words[0]!r}")
    elif not words and not ("text" in kind or "path" in kind):
        errors.append(f"{label}: type cell {cells[1]!r} should say 'text' or 'path'")
    if not cells[3]:
        errors.append(f"{label}: no meaning")
    errors += check_default(label, param, default_cell)
    return errors


def check_default(label: str, param: click.Parameter, cell: str) -> list[str]:
    if param.required:
        return [] if "required" in cell else [f"{label}: default cell should say 'required'"]
    default = param.default
    if default is UNSET:
        # click 8.3 and later mark a parameter without a default this way
        default = False if isinstance(param, click.Option) and param.is_flag else None
    if isinstance(param, click.Option) and param.is_flag and param.secondary_opts:
        want = param.opts[0] if default else param.secondary_opts[0]
        return [] if f"`{want}`" in cell else [f"{label}: default should be `{want}`"]
    if isinstance(param, click.Option) and param.is_flag:
        want = "on" if default else "off"
        return [] if cell.strip("* ") == want else [f"{label}: default should be '{want}'"]
    if default is None or default == "":
        return [] if cell else [f"{label}: empty default cell"]
    if isinstance(default, bool):
        return [] if str(default).lower() in cell.lower() else [f"{label}: default {default}"]
    if isinstance(default, (int, float)):
        numbers = [float(x) for x in re.findall(r"-?\d+(?:\.\d+)?(?:e-?\d+)?", cell.replace("_", ""))]
        ok = any(abs(x - default) <= 1e-12 * max(1.0, abs(default)) for x in numbers)
        return [] if ok else [f"{label}: default cell {cell!r} should give {default!r}"]
    return [] if f"`{default}`" in cell else [f"{label}: default cell {cell!r} should give `{default}`"]


# ---------------------------------------------------------------- shell blocks


def shell_blocks(text: str) -> list[tuple[int, str]]:
    """(line number, body) of every ```bash block."""
    out, lines, i = [], text.splitlines(), 0
    while i < len(lines):
        if re.match(r"^```(bash|sh|shell|console)\s*$", lines[i]):
            start, body = i + 1, []
            i += 1
            while i < len(lines) and not lines[i].startswith("```"):
                body.append(lines[i])
                i += 1
            out.append((start + 1, "\n".join(body)))
        i += 1
    return out


def command_lines(body: str) -> list[str]:
    """`popmix ...` lines of a block, continuations joined, comments dropped."""
    out, lines, i = [], body.splitlines(), 0
    while i < len(lines):
        text = lines[i].strip()
        while text.endswith("\\") and i + 1 < len(lines):
            i += 1
            text = text[:-1] + " " + lines[i].strip()
        text = re.sub(r"^\$\s+", "", text)
        if text.startswith("popmix "):
            out.append(re.split(r"\s(?:#|\||>|2>|&&|;)\s?", text)[0].strip())
        i += 1
    return out


def flags(command: click.Command) -> set[str]:
    return {o for p in command.params for o in (*p.opts, *p.secondary_opts)} | {"--help"}


def check_command_line(text: str) -> str | None:
    words = shlex.split(text)[1:]
    while words and words[0].startswith("-"):
        flag = words.pop(0).split("=")[0]
        if flag not in flags(cli):
            return f"unknown global option {flag}"
        if flag in ("--threads", "-t") and words:
            words.pop(0)
    if not words or words[0] in ("--help", "--version"):
        return None
    group: click.Command = cli
    while isinstance(group, click.Group):
        name = words.pop(0)
        if name not in group.commands:
            return f"no command {name!r}"
        group = group.commands[name]
    known = flags(group)
    for word in words:
        if word.startswith("-") and not re.fullmatch(r"-\d.*", word) and word.split("=")[0] not in known:
            return f"no option {word}"
    return None


def skill_files() -> list[Path]:
    return [SKILL / "SKILL.md", *sorted(COMMANDS.glob("*.md"))]


# ---------------------------------------------------------------- the tests


def test_every_command_has_exactly_one_file():
    files = {p.stem for p in COMMANDS.glob("*.md")}
    assert sorted(set(ALL) - files) == [], "commands without a doc file"
    assert sorted(files - set(ALL)) == [], "doc files without a command"


@pytest.mark.parametrize("stem", sorted(ALL))
def test_file_names_its_command(stem):
    assert doc(stem).splitlines()[0] == f"# popmix {' '.join(PATHS[stem])}"


@pytest.mark.parametrize("stem", sorted(ALL))
def test_every_parameter_is_documented(stem):
    command, rows = ALL[stem], table_rows(doc(stem))
    errors = []
    for param in command.params:
        cells = row_for(param, rows)
        if cells is None:
            errors.append(f"{stem}: no table row for {'/'.join(sorted(spellings(param)))}")
        else:
            errors += check_row(stem, param, cells)
    assert errors == []


@pytest.mark.parametrize("stem", sorted(ALL))
def test_no_row_for_an_option_that_does_not_exist(stem):
    known = {s for p in ALL[stem].params for s in spellings(p)}
    rows = table_rows(doc(stem).split("## Arguments and options", 1)[-1].split("\n## ", 1)[0])
    stale = [t for cells in rows for t in re.findall(r"`([^`]+)`", cells[0]) if t not in known]
    assert stale == []


@pytest.mark.parametrize("stem", sorted(ALL))
def test_file_has_the_sections(stem):
    text = doc(stem)
    for section in ("## Synopsis", "## Arguments and options", "## Inputs and outputs", "## Examples",
                    "## Reading the result", "## Caveats", "## Reference"):
        assert f"\n{section}" in text, f"{stem}.md has no '{section}' section"


def test_global_options_are_documented():
    rows = table_rows((SKILL / "SKILL.md").read_text(encoding="utf-8"))
    errors = []
    for param in cli.params:
        cells = row_for(param, rows)
        if cells is None:
            errors.append(f"SKILL.md: no row for {'/'.join(sorted(spellings(param)))}")
        else:
            errors += check_row("SKILL.md", param, cells)
    assert errors == []


def test_skill_md_frontmatter_and_index():
    text = (SKILL / "SKILL.md").read_text(encoding="utf-8")
    front = re.match(r"^---\n(.*?)\n---\n", text, re.S)
    assert front, "SKILL.md needs YAML frontmatter"
    assert re.search(r"^name: popmix$", front.group(1), re.M)
    description = re.search(r"^description: (.+)$", front.group(1), re.M)
    assert description and "Use when" in description.group(1)
    for stem in ALL:
        assert f"(commands/{stem}.md)" in text, f"SKILL.md does not link commands/{stem}.md"


def test_skill_states_the_popmix_version():
    text = (SKILL / "SKILL.md").read_text(encoding="utf-8")
    stated = re.search(r"^This skill documents popmix (\S+)\.$", text, re.M)
    assert stated, "SKILL.md must say 'This skill documents popmix <version>.'"
    pyproject = re.search(r'^version = "([^"]+)"', (ROOT / "pyproject.toml").read_text(encoding="utf-8"), re.M).group(1)
    assert stated.group(1) == pm.__version__ == pyproject


@pytest.mark.parametrize("path", skill_files(), ids=lambda p: p.name)
def test_shell_blocks_are_marked_and_use_real_flags(path):
    errors = []
    for line, body in shell_blocks(path.read_text(encoding="utf-8")):
        first = body.splitlines()[0] if body else ""
        if not (first == "# runs" or first.startswith("# not run: ")):
            errors.append(f"{path.name}:{line}: first line must be '# runs' or '# not run: <reason>'")
        for text in command_lines(body):
            problem = check_command_line(text)
            if problem:
                errors.append(f"{path.name}:{line}: {text}: {problem}")
    assert errors == []


def test_every_command_has_an_example_that_runs():
    without = [s for s in ALL if not any(b.startswith("# runs") for _, b in shell_blocks(doc(s)))]
    assert without == []


def example_data(directory: Path) -> None:
    """The data SKILL.md's examples use: family, family6 and the small fixtures."""
    for ext in ("geno", "snp", "ind"):
        shutil.copy(f"{FAMILY}.{ext}", directory / f"family.{ext}")
    for ext in ("geno", "snp"):
        shutil.copy(f"{FAMILY}.{ext}", directory / f"family6.{ext}")
    # French and Yoruba each split into three groups by sample order, so that
    # qpAdm and qpWave have enough populations (a toy, not a real labeling)
    seen: dict[str, int] = {}
    out = []
    for line in Path(f"{FAMILY}.ind").read_text(encoding="utf-8").splitlines():
        sample, sex, pop = line.split()
        if pop in ("French", "Yoruba"):
            seen[pop] = seen.get(pop, 0) + 1
            pop = f"{pop}_{seen[pop] % 3 + 1}"
        out.append(f"{sample} {sex} {pop}")
    (directory / "family6.ind").write_text("\n".join(out) + "\n", encoding="utf-8")
    for path in FIXTURES.iterdir():
        if path.is_file():
            shutil.copy(path, directory / path.name)


def run_env(home: Path) -> dict[str, str]:
    env = dict(os.environ)
    env["PATH"] = f"{Path(sys.executable).parent}{os.pathsep}{env.get('PATH', '')}"
    env["HOME"] = str(home)  # never the real ~/.claude
    env["USERPROFILE"] = str(home)  # Path.home() on Windows
    env["MPLBACKEND"] = "Agg"
    return env


# by full path: on Windows a bare "bash" is System32's WSL launcher, found
# before PATH; the examples need a POSIX bash (Git Bash there)
BASH = shutil.which("bash")
POSIX_BASH = BASH is not None and "system32" not in BASH.lower()


@pytest.mark.skipif(not POSIX_BASH, reason="no POSIX bash")
@pytest.mark.skipif(not Path(f"{FAMILY}.geno").exists(), reason="example data not available")
@pytest.mark.parametrize("path", skill_files(), ids=lambda p: p.name)
def test_examples_run(path, tmp_path):
    blocks = [(n, b) for n, b in shell_blocks(path.read_text(encoding="utf-8")) if b.startswith("# runs")]
    if not blocks:
        pytest.skip("no runnable examples")
    work = tmp_path / "work"
    work.mkdir()
    example_data(work)
    for line, body in blocks:
        r = subprocess.run([BASH, "-e", "-o", "pipefail", "-c", body], cwd=work, env=run_env(tmp_path),
                           capture_output=True, text=True, timeout=300)
        assert r.returncode == 0, f"{path.name}:{line} failed:\n{body}\n--- stdout\n{r.stdout}\n--- stderr\n{r.stderr}"
    if path.stem != "install-skill":
        assert not (tmp_path / ".claude").exists()


# ---------------------------------------------------------------- shipping (N4)


def test_skill_ships_in_the_package():
    shipped = resources.files("popmix") / "skill"
    assert (shipped / "SKILL.md").is_file()
    assert sorted(p.name for p in (shipped / "commands").iterdir() if p.name.endswith(".md")) == \
        sorted(f"{s}.md" for s in ALL)


def tree(root: Path) -> dict[str, bytes]:
    return {str(p.relative_to(root)): p.read_bytes() for p in sorted(root.rglob("*")) if p.is_file()}


def test_install_skill_copies_it(tmp_path):
    dest = tmp_path / "skills" / "popmix"
    r = subprocess.run(["popmix", "install-skill", "--dest", str(dest)], env=run_env(tmp_path),
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert tree(dest) == {k: v for k, v in tree(SKILL).items() if "__pycache__" not in k}
    assert pm.__version__ in r.stdout


def test_install_skill_defaults_to_claude_skills_in_home(tmp_path):
    r = subprocess.run(["popmix", "install-skill"], env=run_env(tmp_path), capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert (tmp_path / ".claude" / "skills" / "popmix" / "SKILL.md").is_file()


def test_install_skill_replaces_an_old_copy_only_with_force(tmp_path):
    dest = tmp_path / "popmix"
    dest.mkdir()
    (dest / "SKILL.md").write_text("---\nname: popmix\n---\nold\n", encoding="utf-8")
    (dest / "stale.md").write_text("stale", encoding="utf-8")
    r = subprocess.run(["popmix", "install-skill", "--dest", str(dest)], env=run_env(tmp_path),
                       capture_output=True, text=True)
    assert r.returncode == 1 and "--force" in r.stderr
    assert (dest / "SKILL.md").read_text(encoding="utf-8").endswith("old\n")
    r = subprocess.run(["popmix", "install-skill", "--dest", str(dest), "--force"], env=run_env(tmp_path),
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert not (dest / "stale.md").exists()
    assert (dest / "SKILL.md").read_bytes() == (SKILL / "SKILL.md").read_bytes()


def test_install_skill_refuses_a_file(tmp_path):
    dest = tmp_path / "popmix"
    dest.write_text("not a directory", encoding="utf-8")
    r = subprocess.run(["popmix", "install-skill", "--dest", str(dest), "--force"], env=run_env(tmp_path),
                       capture_output=True, text=True)
    assert r.returncode == 1
    assert dest.read_text(encoding="utf-8") == "not a directory"


def test_install_skill_never_replaces_something_else(tmp_path):
    dest = tmp_path / "notes"
    dest.mkdir()
    (dest / "mine.txt").write_text("keep", encoding="utf-8")
    r = subprocess.run(["popmix", "install-skill", "--dest", str(dest), "--force"], env=run_env(tmp_path),
                       capture_output=True, text=True)
    assert r.returncode == 1 and "not a popmix skill" in r.stderr
    assert (dest / "mine.txt").read_text(encoding="utf-8") == "keep"
