"""How the CLI reports errors: expected ones in one line, bugs in full."""

from __future__ import annotations

import os
import signal
import sys
from typing import NoReturn

import click

import popmix as pm

# Failures of the input (bad names or values, missing or malformed files) or
# of the computation (a singular system): reported in one line. Anything else
# is a bug and keeps its traceback.
EXPECTED_ERRORS = (ValueError, OSError, RuntimeError)


def fail(error: Exception) -> NoReturn:
    """Report an expected error and exit 1; with --debug, raise it.

    A popmix InternalError is a bug, not an expected error: it is raised,
    with its traceback, as other bugs are.
    """
    if isinstance(error, pm.InternalError):
        raise error
    if isinstance(error, BrokenPipeError):

        # `popmix ... | head`: the reader closed the pipe. Exit as a process
        # killed by SIGPIPE, without a message, with stdout on /dev/null so
        # the interpreter's final flush cannot fail again.
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        sys.exit(128 + signal.SIGPIPE)
    ctx = click.get_current_context(silent=True)
    if ctx is not None and (ctx.find_root().obj or {}).get("debug"):
        raise error
    click.echo(f"Error: {error}", err=True)
    sys.exit(1)
