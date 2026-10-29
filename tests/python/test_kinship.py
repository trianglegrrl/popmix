"""Kinship API and the kinship command, on the 54-sample DATES data.

The estimates are checked against PLINK 2 (KING-robust) and READv2 in
tests/golden/test_king_plink2.py and tests/golden/test_read_readv2.py.
"""

import os
import subprocess
import sys

import pandas as pd
import pytest

import popmix as pm

SMALL = "tests/golden/dates/family_chr20_22"

pytestmark = pytest.mark.skipif(not os.path.exists(f"{SMALL}.geno"), reason="test data not available")


def test_king_on_all_populations():
    df = pm.king(SMALL)
    assert len(df) == 54 * 53 // 2
    assert list(df.columns) == ["id1", "id2", "nsnp", "hethet", "ibs0", "het1_hom2", "het2_hom1", "kinship"]


def test_read_normalizes_by_the_median_pair():
    df = pm.read_kinship(SMALL, ["French", "Yoruba"])
    assert df["normalized_p0"].median() == pytest.approx(1.0)
    assert df.attrs["norm_value"] == pytest.approx(df["p0"].median())


def test_read_with_a_given_normalization():
    df = pm.read_kinship(SMALL, ["French"], norm_value=0.25)
    assert df.attrs["norm_value"] == 0.25
    assert (df["normalized_p0"] == df["p0"] / 0.25).all()


def test_read_rejects_unknown_options():
    with pytest.raises(ValueError, match="norm must be"):
        pm.read_kinship(SMALL, ["French"], norm="mode")
    with pytest.raises(ValueError, match="unknown compat"):
        pm.read_kinship(SMALL, ["French"], compat="readv1")


@pytest.mark.parametrize("method", ["king", "read"])
def test_kinship_command_matches_the_api(tmp_path, method):
    out = tmp_path / "kin.tsv"
    subprocess.run(["popmix", "kinship", "--input", SMALL, "--populations", "French,Yoruba", "--method", method,
                    "--output", str(out)], check=True, capture_output=True, text=True)
    got = pd.read_csv(out, sep="\t")
    api = (pm.king if method == "king" else pm.read_kinship)(SMALL, ["French", "Yoruba"])
    assert got["id1"].tolist() == api["id1"].tolist()
    assert got["kinship"].to_numpy() == pytest.approx(api["kinship"].to_numpy(), rel=1e-11)


def test_columns_are_written_as_the_table_writer_writes_them(tmp_path):
    """The KING command writes its table without pandas, in the same format."""
    import numpy as np

    from popmix.cli import write_columns, write_table

    cols = {"id1": ["a", "b", "c"], "nsnp": np.array([3, 0, 12], dtype=np.uint64),
            "kinship": np.array([0.1234567890123456, -np.inf, np.nan]), "x": np.array([1e-20, 2.5, -0.0])}
    write_table(pd.DataFrame(cols), "tsv", str(tmp_path / "pandas.tsv"))
    write_columns(cols, str(tmp_path / "columns.tsv"))
    assert (tmp_path / "columns.tsv").read_text(encoding="utf-8") == (tmp_path / "pandas.tsv").read_text(encoding="utf-8")


def test_the_king_command_does_not_load_pandas():
    code = ("import sys\nfrom click.testing import CliRunner\nfrom popmix.cli import cli\n"
            f"r = CliRunner().invoke(cli, ['kinship', '-i', '{SMALL}', '--populations', 'French', '--method', 'king'])\n"
            "assert r.exit_code == 0, r.output\nprint('pandas' in sys.modules)")
    out = subprocess.run([sys.executable, "-c", code], check=True, capture_output=True, text=True).stdout
    assert out.strip() == "False"


def test_tables_end_their_lines_with_a_newline_only(tmp_path):
    """On every platform: pandas would end lines with os.linesep, which a text
    file on Windows turns into a blank line after every row."""
    from popmix.cli import write_table

    write_table(pd.DataFrame({"a": [1, 2]}), "tsv", str(tmp_path / "t.tsv"))
    assert (tmp_path / "t.tsv").read_bytes() == b"a\n1\n2\n"
