"""The KING table written in Rust is the text write_columns writes (TASKS.md F14).

``popmix kinship --method king`` builds its TSV in Rust (``king_tsv``), so
that the command imports neither numpy nor pandas and loops over no pair in
Python; the text must be the same as writing ``king_columns`` with
``write_columns``, numbers at ``%.12g``.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from popmix.cli._common import write_columns
from popmix.kinship import king_columns, king_tsv

ROOT = Path(__file__).parents[2]
SMALL = str(ROOT / "tests" / "golden" / "dates" / "family_chr20_22")


def test_the_rust_table_is_write_columns_text(tmp_path):
    out = tmp_path / "cols.tsv"
    write_columns(king_columns(SMALL, None), str(out))
    assert king_tsv(SMALL) == out.read_text(encoding="utf-8")


def test_a_minus_infinite_kinship_is_written_as_python_writes_it(tmp_path):
    # a pseudo-haploid sample (no heterozygous call) has kinship -inf
    prefix = str(tmp_path / "p")
    subprocess.run(["popmix", "convertf", "-i", SMALL, "-o", prefix, "-of", "eigenstrat"], check=True,
                   capture_output=True)
    lines = Path(f"{prefix}.geno").read_text(encoding="utf-8").splitlines()
    Path(f"{prefix}.geno").write_text("".join(("2" if x[0] == "1" else x[0]) + x[1:] + "\n" for x in lines), encoding="utf-8")
    out = tmp_path / "cols.tsv"
    write_columns(king_columns(prefix, ["Simulation", "French"]), str(out))
    text = king_tsv(prefix, ["Simulation", "French"])
    assert "-inf" in text
    assert text == out.read_text(encoding="utf-8")


def test_the_command_imports_no_numpy(tmp_path):
    code = ("import sys\nfrom popmix.cli import cli\n"
            f"sys.argv = ['popmix', 'kinship', '-i', {SMALL!r}, '-o', {str(tmp_path / 'k.tsv')!r}]\n"
            "try:\n    cli()\nexcept SystemExit:\n    pass\nprint('numpy' in sys.modules)")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert r.stdout.strip().splitlines()[-1] == "False"
