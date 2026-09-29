"""The stages of qdex.cli.main declare their inputs; each must be exported by an earlier stage."""
import ast
import inspect
from pathlib import Path

import qdex.cli as cli

STAGES = [cli._prepare_run, cli._read_geometry_and_basis, cli._read_molecular_orbitals,
          cli._quasiparticle_correction, cli._qp_levels_and_kernel, cli._ip_ea_and_energy_axis,
          cli._orbital_populations, cli._active_space_and_soc, cli._cubes_and_fuzzy,
          cli._transition_dipoles, cli._solve_excitons]


def _exports(func):
    """Names in the stage's final ``return _export(locals(), (...))``."""
    tree = ast.parse(inspect.getsource(func))
    names = set()
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and getattr(node.func, "id", None) == "_export"
                and isinstance(node.args[1], ast.Tuple)):
            names |= {elt.value for elt in node.args[1].elts}
    return names


def test_main_runs_the_listed_stages():
    source = inspect.getsource(cli.main)
    for stage in STAGES:
        assert stage.__name__ in source


def test_every_stage_input_is_exported_earlier():
    available = {"parser", "config_path"}
    for stage in STAGES:
        params = [p.name for p in inspect.signature(stage).parameters.values() if p.kind is p.KEYWORD_ONLY]
        missing = sorted(set(params) - available)
        assert not missing, f"{stage.__name__} needs {missing}, which no earlier stage exports"
        available |= _exports(stage)


def test_stage_variables_are_declared():
    """A name a stage reads but neither defines nor declares would fail at run time."""
    import builtins
    module_names = set(dir(cli)) | set(dir(builtins))
    for stage in STAGES:
        tree = ast.parse(inspect.getsource(stage))
        func = tree.body[0]
        declared = {a.arg for a in func.args.args + func.args.kwonlyargs}
        stored = {n.id for n in ast.walk(func) if isinstance(n, ast.Name) and isinstance(n.ctx, (ast.Store, ast.Del))}
        stored |= {n.name for n in ast.walk(func) if isinstance(n, (ast.FunctionDef, ast.ClassDef))}
        stored |= {(a.asname or a.name).split(".")[0] for n in ast.walk(func)
                   if isinstance(n, (ast.Import, ast.ImportFrom)) for a in n.names}
        stored |= {n.name for n in ast.walk(func) if isinstance(n, ast.ExceptHandler) and n.name}
        stored |= {a.arg for n in ast.walk(func) if isinstance(n, (ast.FunctionDef, ast.Lambda))
                   for a in n.args.args + n.args.kwonlyargs}
        loaded = {n.id for n in ast.walk(func) if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)}
        unknown = sorted(loaded - declared - stored - module_names)
        assert not unknown, f"{stage.__name__} reads undeclared names {unknown}"
