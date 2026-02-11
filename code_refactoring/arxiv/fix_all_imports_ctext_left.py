#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check and fix ALL Python scripts in a directory against felics_renames.json.

- Reads felics_renames.json (produced by the defs pipeline).
- Walks every .py file under the target directory (recursively).
- For each rename (old -> new), checks if the file still uses the old name.
- If so, rewrites imports/usages to the new name.
- Reports what was changed per file.

Import behavior:
- `from module import Name`:
    * `module` is left unchanged.
    * `Name` (and any alias) is renamed if in the rename_map.
- `import Name`:
    * `Name` (and any alias) is renamed if in the rename_map.
"""

# Standard libraries
import json
import pathlib
import sys

# Third party libraries
import libcst as cst
from libcst import metadata


class UsageFixer(cst.CSTTransformer):
    def __init__(
        self,rename_map,):
        self.rename_map = rename_map
        self.changes = []

    # --- from module import Name ---
    # Only rename the imported names, NOT the module part.
    def leave__import_from(
        self,orig,updated,):
        # Don't touch "from module import *"
        if isinstance(
            updated.names,cst.ImportStar,):
            return updated

        new_aliases = []
        for alias in updated.names:
            if isinstance(
                alias,cst.ImportAlias,):
                # from module import Name as Alias
                if alias.asname:
                    old = alias.asname.name.value
                    new = self.rename_map.get(
                        old,old,)
                    if new != old:
                        self.changes.append((old, new))
                        alias = alias.with_changes(
                            asname=cst.AsName(name=cst.Name(new))
                        )
                # from module import Name
                else:
                    old = alias.name.value
                    new = self.rename_map.get(
                        old,old,)
                    if new != old:
                        self.changes.append((old, new))
                        alias = alias.with_changes(name=cst.Name(new))
            new_aliases.append(alias)

        # module part (updated.module) is left untouched
        return updated.with_changes(names=tuple(new_aliases))

    # --- import Name ---
    # Rename either the imported name or its alias.
    def leave__import(
        self,orig,updated,):
        new_names = []
        for alias in updated.names:
            if isinstance(
                alias,cst.ImportAlias,):
                # import Name as Alias
                if alias.asname:
                    old = alias.asname.name.value
                    new = self.rename_map.get(
                        old,old,)
                    if new != old:
                        self.changes.append((old, new))
                        alias = alias.with_changes(
                            asname=cst.AsName(name=cst.Name(new))
                        )
                # import Name
                else:
                    old = alias.name.value
                    new = self.rename_map.get(
                        old,old,)
                    if new != old:
                        self.changes.append((old, new))
                        alias = alias.with_changes(name=cst.Name(new))

            new_names.append(alias)
        return updated.with_changes(names=tuple(new_names))

    # --- plain names in code (variables, functions, etc.) ---
    def leave__name(
        self,orig,updated,):
        old = updated.value
        new = self.rename_map.get(
            old,old,)
        if new != old:
            self.changes.append((old, new))
            return cst.Name(new)
        return updated

    # --- attributes like obj.old_name ---
    def leave__attribute(
        self,orig,updated,):
        old_attr = updated.attr.value
        new_attr = self.rename_map.get(
            old_attr,old_attr,)
        if new_attr != old_attr:
            self.changes.append((old_attr, new_attr))
            return updated.with_changes(attr=cst.Name(new_attr))
        return updated


def main():
    if len(sys.argv) != 3:
        print(
            "Usage: felics_check_and_fix_all.py <felics_renames.json> <target_dir>"
        )
        sys.exit(1)

    log_path = pathlib.Path(sys.argv[1]).resolve()
    target_dir = pathlib.Path(sys.argv[2]).resolve()

    if not log_path.exists():
        print(f"Rename log not found: {log_path}")
        sys.exit(1)
    if not target_dir.exists():
        print(f"Target directory not found: {target_dir}")
        sys.exit(1)

    entries = json.loads(log_path.read_text(encoding="utf-8"))
    # entries are [file_path, old, new]
    rename_map = {old: new for _, old, new in entries}

    print(f"Using rename map with {len(rename_map)} entries.")
    print(f"Scanning directory: {target_dir}")

    for pyfile in target_dir.rglob("*.py"):
        try:
            source = pyfile.read_text(encoding="utf-8")
            module = cst.parse_module(source)
        except Exception:
            continue

        wrapper = metadata.MetadataWrapper(module)
        fixer = UsageFixer(rename_map)
        new_tree = wrapper.visit(fixer)

        if fixer.changes:
            print(f"\nApplied {len(fixer.changes)} changes in {pyfile}:")
            seen = set()
            for old, new in fixer.changes:
                if (old, new) not in seen:
                    seen.add((old, new))
                    print(f"  {old} -> {new}")
            pyfile.write_text(
                new_tree.code,encoding="utf-8",)
        # else: silent if no changes


if __name__ == "__main__":
    main()

# usage
# python code_refactoring/felics_check_and_fix_all.py code_refactoring/felics_renames.json .
