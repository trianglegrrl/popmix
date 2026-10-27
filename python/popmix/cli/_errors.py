"""How the CLI reports errors: expected ones in one line, bugs in full."""

from __future__ import annotations

import os
import signal
import sys
import warnings
from typing import NoReturn

import click

import popmix as pm

# Failures of the input (bad names or values, missing or malformed files) or
# of the computation (a singular system): reported in one line. Anything else
# is a bug and keeps its traceback.
EXPECTED_ERRORS = (ValueError, OSError, RuntimeError)


def print_popmix_warnings() -> None:
    """Show each PopmixWarning as one ``Warning: ...`` line on stderr, without
    Python's file, line and source; other warnings as Python shows them."""
    show = warnings.showwarning

    def showwarning(message, category, filename, lineno, file=None, line=None):
        if issubclass(category, pm.PopmixWarning):
            click.echo(f"Warning: {message}", err=True)
        else:
            show(message, category, filename, lineno, file, line)

    warnings.showwarning = showwarning


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
