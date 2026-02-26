#!/usr/bin/env python3
"""
FELiCS two-pass style rewriter

Pipeline A ("defs"):
  - Lift imports to top-level
  - PascalCase classes
  - snake_case functions
  - Split args/params one per line + trailing commas for multi-arg calls/defs
  - Remove duplicate top-level imports after lifting
  - Log renames to felics_renames.json

Pipeline B ("usages"):
  - Read felics_renames.json
  - Apply renames consistently across imports/usages
"""

import json
import pathlib
import re
import sys
from typing import Dict, List, Tuple

import libcst as cst
from libcst import metadata

SCRIPT_DIR = "/Users/marina/Documents/FELiCS-main/code_refactoring/"

# ----------------------------
# Helpers
# ----------------------------


def to_snake(name: str) -> str:
    prefix = ""
    while name.startswith("_"):
        prefix += "_"
        name = name[1:]
    s1 = re.sub(r"(.)([A-Z][a-z]+)", r"\1_\2", name)
    s2 = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", s1)
    snake = s2.lower()
    return prefix + snake


def to_pascal(name: str) -> str:
    prefix = ""
    while name.startswith("_"):
        prefix += "_"
        name = name[1:]
    if not name:
        return prefix
    # timeStepping -> TimeStepping
    if re.match(r"[a-z].*[A-Z]", name):
        return prefix + name[0].upper() + name[1:]
    # Already Camel/Pascal (has another capital, starts upper) -> keep
    if re.match(r"[A-Z].*[A-Z]", name):
        return prefix + name
    # else just uppercase first char
    return prefix + name[0].upper() + name[1:]


def _one_per_line_ws() -> cst.ParenthesizedWhitespace:
    # newline then indent one level inside parentheses
    return cst.ParenthesizedWhitespace(
        first_line=cst.TrailingWhitespace(
            whitespace=cst.SimpleWhitespace(""),
            newline=cst.Newline(),
        ),
        indent=True,
    )


def _split_call_args_one_per_line(call: cst.Call) -> cst.Call:
    if len(call.args) <= 1:
        return call

    ws = _one_per_line_ws()

    # 1) newline right after "(" for the first argument
    # LibCST supports this via whitespace_before_args on many versions.
    if hasattr(call, "whitespace_before_args"):
        call = call.with_changes(whitespace_before_args=ws)

    # 2) newline+indent after each comma, and force trailing comma
    new_args = [arg.with_changes(comma=cst.Comma(whitespace_after=ws)) for arg in call.args]
    return call.with_changes(args=new_args)


def _split_parameters_one_per_line(params: cst.Parameters) -> cst.Parameters:
    # Count params including *args/**kwargs (ignore empty/None)
    count = len(params.posonly_params) + len(params.params) + len(params.kwonly_params)
    if isinstance(params.star_arg, cst.Param):
        count += 1
    if isinstance(params.star_kwarg, cst.Param):
        count += 1

    if count <= 1:
        return params

    ws = _one_per_line_ws()

    # newline right after "(" in defs if we can safely modify existing lpar
    # (DO NOT create new lpar/rpar; that can change semantics)
    if hasattr(params, "lpar") and isinstance(params.lpar, tuple) and len(params.lpar) > 0:
        new_lpar0 = params.lpar[0].with_changes(whitespace_after=ws)
        params = params.with_changes(lpar=(new_lpar0,) + params.lpar[1:])

    def add_commas(seq):
        return [item.with_changes(comma=cst.Comma(whitespace_after=ws)) for item in seq] if seq else seq

    new_star_arg = params.star_arg
    if isinstance(new_star_arg, cst.Param):
        new_star_arg = new_star_arg.with_changes(comma=cst.Comma(whitespace_after=ws))

    new_star_kwarg = params.star_kwarg
    if isinstance(new_star_kwarg, cst.Param):
        new_star_kwarg = new_star_kwarg.with_changes(comma=cst.Comma(whitespace_after=ws))

    return params.with_changes(
        posonly_params=add_commas(params.posonly_params),
        params=add_commas(params.params),
        star_arg=new_star_arg,
        kwonly_params=add_commas(params.kwonly_params),
        star_kwarg=new_star_kwarg,
    )


def _dedupe_top_level_imports(module: cst.Module) -> cst.Module:
    """
    Remove exact-duplicate import lines at the top of the module.
    This only de-dupes *top-level contiguous import statements* (common after lifting).
    """
    body = list(module.body)
    if not body:
        return module

    def is_import_stmt(stmt: cst.CSTNode) -> bool:
        if not isinstance(stmt, cst.SimpleStatementLine):
            return False
        return any(isinstance(elem, (cst.Import, cst.ImportFrom)) for elem in stmt.body)

    i = 0
    # keep leading module docstring untouched
    if (
        len(body) >= 1
        and isinstance(body[0], cst.SimpleStatementLine)
        and len(body[0].body) == 1
        and isinstance(body[0].body[0], cst.Expr)
        and isinstance(body[0].body[0].value, cst.SimpleString)
    ):
        i = 1

    # De-dupe only contiguous top-level imports starting at i
    seen: set[str] = set()
    new_body = body[:i]
    while i < len(body) and is_import_stmt(body[i]):
        stmt = body[i]
        key = module.code_for_node(stmt).strip()
        if key not in seen:
            seen.add(key)
            new_body.append(stmt)
        i += 1

    new_body.extend(body[i:])
    return module.with_changes(body=tuple(new_body))


# ----------------------------
# Import lifting
# ----------------------------
class ImportLifter(cst.CSTTransformer):
    METADATA_DEPENDENCIES = (metadata.ParentNodeProvider,)

    def __init__(self) -> None:
        self.lifted_statements: List[cst.SimpleStatementLine] = []

    def leave_SimpleStatementLine(self, orig, updated):
        has_import = any(
            isinstance(elem, (cst.Import, cst.ImportFrom))
            for elem in updated.body
        )
        if not has_import:
            return updated
        parent = self.get_metadata(metadata.ParentNodeProvider, orig)
        if not isinstance(parent, cst.Module):
            self.lifted_statements.append(updated)
            return cst.RemoveFromParent()
        return updated


# ----------------------------
# Pipeline A: Definition renamer
# ----------------------------
class DefinitionRenamer(cst.CSTTransformer):
    def __init__(
        self,
        file_path: pathlib.Path,
        renames_log: List[Tuple[str, str, str]],
    ) -> None:
        self.file_path = str(file_path)
        self.renames_log = renames_log

    def leave_ClassDef(self, orig, updated):
        old = updated.name.value
        new = to_pascal(old)
        if new != old:
            self.renames_log.append((self.file_path, old, new))
        return updated.with_changes(name=cst.Name(new))

    def leave_FunctionDef(self, orig, updated):
        old = updated.name.value
        new = to_snake(old)
        if new != old:
            self.renames_log.append((self.file_path, old, new))
        return updated.with_changes(name=cst.Name(new))

    def leave_Call(self, orig, updated):
        # Split call args one per line + trailing commas
        return _split_call_args_one_per_line(updated)

    def leave_Parameters(self, orig, updated):
        # Split parameters one per line + trailing commas
        return _split_parameters_one_per_line(updated)


# ----------------------------
# Pipeline B: Usage renamer
# ----------------------------
class UsageRenamer(cst.CSTTransformer):
    def __init__(self, rename_map: Dict[str, str]) -> None:
        self.rename_map = rename_map

    def leave_ImportFrom(self, orig, updated):
        if isinstance(updated.names, cst.ImportStar):
            return updated
        new_aliases = []
        for alias in updated.names:
            if isinstance(alias, cst.ImportAlias):
                # from module import Name as Alias
                if alias.asname:
                    old = alias.asname.name.value
                    new = self.rename_map.get(old, old)
                    if new != old:
                        alias = alias.with_changes(
                            asname=cst.AsName(name=cst.Name(new))
                        )
                # from module import Name
                else:
                    old = alias.name.value
                    new = self.rename_map.get(old, old)
                    if new != old:
                        alias = alias.with_changes(name=cst.Name(new))
            new_aliases.append(alias)
        return updated.with_changes(names=tuple(new_aliases))

    def leave_Import(self, orig, updated):
        new_names = []
        for alias in updated.names:
            if isinstance(alias, cst.ImportAlias) and alias.asname:
                old = alias.asname.name.value
                new = self.rename_map.get(old, old)
                if new != old:
                    alias = alias.with_changes(
                        asname=cst.AsName(name=cst.Name(new))
                    )
            new_names.append(alias)
        return updated.with_changes(names=tuple(new_names))

    def leave_Name(self, orig, updated):
        old = updated.value
        new = self.rename_map.get(old, old)
        if new != old:
            return cst.Name(new)
        return updated

    def leave_Attribute(self, orig, updated):
        old_attr = updated.attr.value
        new_attr = self.rename_map.get(old_attr, old_attr)
        if new_attr != old_attr:
            return updated.with_changes(attr=cst.Name(new_attr))
        return updated

    def leave_Call(self, orig, updated):
        # Split call args one per line + trailing commas
        return _split_call_args_one_per_line(updated)


# ----------------------------
# Orchestration
# ----------------------------
def pipeline_a(root: pathlib.Path) -> None:
    renames_log: List[Tuple[str, str, str]] = []
    for pyfile in root.rglob("*.py"):
        try:
            source = pyfile.read_text(encoding="utf-8")
            module = cst.parse_module(source)
        except Exception:
            continue

        wrapper = metadata.MetadataWrapper(module)
        import_lifter = ImportLifter()
        lifted_tree = wrapper.visit(import_lifter)

        if import_lifter.lifted_statements:
            existing_body = list(lifted_tree.body)
            new_body = list(import_lifter.lifted_statements) + existing_body
            lifted_tree = lifted_tree.with_changes(body=tuple(new_body))

        # Remove exact duplicate top-level import lines after lifting
        lifted_tree = _dedupe_top_level_imports(lifted_tree)

        transformer = DefinitionRenamer(pyfile, renames_log)
        new_tree = lifted_tree.visit(transformer)
        pyfile.write_text(new_tree.code, encoding="utf-8")

    log_path = pathlib.Path(SCRIPT_DIR) / "felics_renames.json"
    log_path.write_text(
        json.dumps(renames_log, indent=2),
        encoding="utf-8",
    )
    print(
        f"Pipeline A complete. Wrote {len(renames_log)} renames to {log_path}"
    )


def pipeline_b(root: pathlib.Path) -> None:
    log_path = pathlib.Path(SCRIPT_DIR) / "felics_renames.json"
    if not log_path.exists():
        print(
            "No renames log found. Run: felics_style_rewriter.py defs <target_dir>"
        )
        return
    try:
        entries: List[Tuple[str, str, str]] = json.loads(
            log_path.read_text(encoding="utf-8")
        )
    except Exception as e:
        print(f"Failed to read renames log: {e}")
        return

    rename_map: Dict[str, str] = {}
    for _file, old, new in entries:
        if old != new:
            rename_map[old] = new

    for pyfile in root.rglob("*.py"):
        try:
            source = pyfile.read_text(encoding="utf-8")
            module = cst.parse_module(source)
        except Exception:
            continue
        transformer = UsageRenamer(rename_map)
        new_tree = module.visit(transformer)
        pyfile.write_text(new_tree.code, encoding="utf-8")

    print("Pipeline B complete. Applied renames across imports and usages.")


def main() -> None:
    if len(sys.argv) != 3 or sys.argv[1] not in {"defs", "usages"}:
        print("Usage: felics_style_rewriter.py [defs|usages] <target_dir>")
        sys.exit(1)
    mode = sys.argv[1]
    root = pathlib.Path(sys.argv[2]).resolve()
    if mode == "defs":
        pipeline_a(root)
    else:
        pipeline_b(root)


if __name__ == "__main__":
    main()
