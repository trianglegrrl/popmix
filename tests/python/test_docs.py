"""Docstrings and type hints of the public API (TASKS.md J3, CLAUDE.md).

Every public function, class and method has a NumPy-format docstring
(Parameters when it takes any, Returns), and the compiled module has a
type stub, ``popmix/_core.pyi``, that declares every public name of it
with the parameters the module really takes.
"""

from __future__ import annotations

import ast
import inspect
from pathlib import Path

import pytest

import popmix as pm

STUB = Path(pm.__file__).with_name("_core.pyi")
NO_RETURNS = ("__repr__", "__reduce__")


def _members():
    for name in pm.__all__:
        obj = getattr(pm, name)
        if isinstance(obj, type):
            if issubclass(obj, (BaseException, Warning)):
                continue
            yield name, obj
            for attr in dir(obj):
                if not attr.startswith("_") and callable(getattr(obj, attr)):
                    yield f"{name}.{attr}", getattr(obj, attr)
        elif callable(obj):
            yield name, obj


def _params(obj) -> list[str] | None:
    try:
        return [p for p in inspect.signature(obj).parameters if p not in ("self", "cls")]
    except (TypeError, ValueError):
        return None


@pytest.mark.parametrize("name,obj", list(_members()), ids=lambda x: x if isinstance(x, str) else "")
def test_public_names_have_numpy_docstrings(name, obj):
    doc = inspect.getdoc(obj) or ""
    assert doc, "no docstring"
    if isinstance(obj, type):
        return
    if _params(obj):
        assert "Parameters\n----------" in doc
    if not name.endswith(NO_RETURNS):
        assert "Returns\n-------" in doc


def _stub() -> dict[str, ast.AST]:
    tree = ast.parse(STUB.read_text(encoding="utf-8"))
    return {n.name: n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.ClassDef))}


def _stub_params(fn: ast.FunctionDef) -> list[str]:
    a = fn.args
    return [x.arg for x in a.posonlyargs + a.args + a.kwonlyargs if x.arg not in ("self", "cls")]


def _core_public():
    for name in pm.__all__:
        obj = getattr(pm, name)
        if getattr(obj, "__module__", "") == "popmix._core" and not name.startswith("__"):
            yield name, obj


def test_the_compiled_module_has_a_type_stub_and_marker():
    assert STUB.exists()
    assert STUB.with_name("py.typed").exists()


@pytest.mark.parametrize("name,obj", list(_core_public()), ids=lambda x: x if isinstance(x, str) else "")
def test_the_stub_declares_what_the_module_takes(name, obj):
    stub = _stub()
    assert name in stub, f"{name} missing from _core.pyi"
    if isinstance(obj, type):
        methods = {n.name: n for n in stub[name].body if isinstance(n, ast.FunctionDef)}
        for attr in dir(obj):
            if attr.startswith("_") or not callable(getattr(obj, attr)):
                continue
            assert attr in methods, f"{name}.{attr} missing from _core.pyi"
            assert _stub_params(methods[attr]) == _params(getattr(obj, attr)), f"{name}.{attr}"
    else:
        assert _stub_params(stub[name]) == _params(obj)
