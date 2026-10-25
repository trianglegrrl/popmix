"""Every tracked file can be checked out on Windows (the release builds there).

Windows refuses the characters < > : " | ? * and names ending in a dot or
space; DATES writes `name.out:1` and `name:log`.
"""

import re
import shutil
import subprocess

import pytest

FORBIDDEN = re.compile(r'[<>:"|?*\\]|[. ]$')


@pytest.mark.skipif(shutil.which("git") is None, reason="git not available")
def test_tracked_names_are_portable():
    run = subprocess.run(["git", "ls-files", "-z"], capture_output=True, text=True, check=True)
    bad = [path for path in run.stdout.split("\0") if path and any(FORBIDDEN.search(part) for part in path.split("/"))]
    assert bad == []
