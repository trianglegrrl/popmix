"""The exceptions popmix raises.

Every one is a :class:`PopmixError`, and also the built-in exception the
same failure would be in plain Python, so code that catches ``ValueError``
or ``FileNotFoundError`` keeps working::

    PopmixError
    ├── InputError (ValueError)             a bad argument, file or name
    │   └── PopulationNotFoundError         a population not in the data
    ├── MissingFileError (FileNotFoundError)
    ├── FileError (OSError)                 reading or writing a file failed
    ├── ComputationError (RuntimeError)     valid input, failed computation
    └── InternalError (RuntimeError)        a bug in popmix: please report it
"""

from __future__ import annotations


class PopmixError(Exception):
    """Base class of every error popmix raises."""


class InputError(PopmixError, ValueError):
    """The input is wrong: an argument, a malformed file, a name."""


class PopulationNotFoundError(InputError):
    """A population that is not in the data or the f2 blocks."""


class MissingFileError(PopmixError, FileNotFoundError):
    """A file that is not there."""


class FileError(PopmixError, OSError):
    """Reading or writing a file failed."""


class ComputationError(PopmixError, RuntimeError):
    """The computation failed on valid input (a singular system, say)."""


class InternalError(PopmixError, RuntimeError):
    """A bug in popmix, not a problem with the input: please report it."""


__all__ = [
    "PopmixError",
    "InputError",
    "PopulationNotFoundError",
    "MissingFileError",
    "FileError",
    "ComputationError",
    "InternalError",
]
