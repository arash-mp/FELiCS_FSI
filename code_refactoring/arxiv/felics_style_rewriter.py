#!/usr/bin/env python3
"""
FELiCS two-pass style rewriter

Pipeline A ("defs"):
  - Lift imports to top-level
  - PascalCase classes
  - snake_case functions
  - Trailing commas for multi-arg calls/defs
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
        if len(updated.args) > 1:
            new_args = [
                arg.with_changes(comma=cst.Comma()) for arg in updated.args
            ]
            return updated.with_changes(args=new_args)
        return updated

    def leave_Parameters(self, orig, updated):
        def add_commas(seq):
            return (
                [item.with_changes(comma=cst.Comma()) for item in seq]
                if seq
                else seq
            )

        return updated.with_changes(
            posonly_params=add_commas(updated.posonly_params),
            params=add_commas(updated.params),
            kwonly_params=add_commas(updated.kwonly_params),
        )


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
        if len(updated.args) > 1:
            new_args = [
                arg.with_changes(comma=cst.Comma()) for arg in updated.args
            ]
            updated = updated.with_changes(args=new_args)
        return updated


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
