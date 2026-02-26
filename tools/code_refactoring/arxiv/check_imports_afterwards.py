#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 12:02:00 2025


Check and fix a single Python script against felics_renames.json.

- Reads felics_renames.json (produced by the defs pipeline).
- Loads a target .py file.
- For each rename (old -> new), checks if the file still uses the old name.
- If so, rewrites imports/usages to the new name.
- Reports what was changed.
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
            self, rename_map,):
        self.rename_map = rename_map
        self.changes = []

    def leave__import_from(
            self, orig, updated,):
        if isinstance(
                updated.names, cst.ImportStar,):
            return updated
        new_aliases = []
        for alias in updated.names:
            if isinstance(
                    alias, cst.ImportAlias,):
                if alias.asname:
                    old = alias.asname.name.value
                    new = self.rename_map.get(
                        old, old,)
                    if new != old:
                        self.changes.append((old, new))
                        alias = alias.with_changes(
                            asname=cst.AsName(name=cst.Name(new))
                        )
                else:
                    old = alias.name.value
                    new = self.rename_map.get(
                        old, old,)
                    if new != old:
                        self.changes.append((old, new))
                        alias = alias.with_changes(name=cst.Name(new))
            new_aliases.append(alias)
        return updated.with_changes(names=tuple(new_aliases))

    def leave__import(
            self, orig, updated,):
        new_names = []
        for alias in updated.names:
            if (
                isinstance(
                    alias, cst.ImportAlias,)
                and alias.asname
            ):
                old = alias.asname.name.value
                new = self.rename_map.get(
                    old, old,)
                if new != old:
                    self.changes.append((old, new))
                    alias = alias.with_changes(
                        asname=cst.AsName(name=cst.Name(new))
                    )
            new_names.append(alias)
        return updated.with_changes(names=tuple(new_names))

    def leave__name(
            self, orig, updated,):
        old = updated.value
        new = self.rename_map.get(
            old, old,)
        if new != old:
            self.changes.append((old, new))
            return cst.Name(new)
        return updated

    def leave__attribute(
            self, orig, updated,):
        old_attr = updated.attr.value
        new_attr = self.rename_map.get(
            old_attr, old_attr,)
        if new_attr != old_attr:
            self.changes.append((old_attr, new_attr))
            return updated.with_changes(attr=cst.Name(new_attr))
        return updated


def main():
    if len(sys.argv) != 3:
        print(
            "Usage: felics_check_and_fix.py <felics_renames.json> <target_file.py>"
        )
        sys.exit(1)

    log_path = pathlib.Path(sys.argv[1])
    target_file = pathlib.Path(sys.argv[2])

    if not log_path.exists():
        print(f"Rename log not found: {log_path}")
        sys.exit(1)
    if not target_file.exists():
        print(f"Target file not found: {target_file}")
        sys.exit(1)

    entries = json.loads(log_path.read_text(encoding="utf-8"))
    rename_map = {old: new for _, old, new in entries}

    source = target_file.read_text(encoding="utf-8")
    module = cst.parse_module(source)
    wrapper = metadata.MetadataWrapper(module)

    fixer = UsageFixer(rename_map)
    new_tree = wrapper.visit(fixer)

    if fixer.changes:
        print(f"Applied {len(fixer.changes)} changes in {target_file}:")
        for old, new in fixer.changes:
            print(f"  {old} -> {new}")
        target_file.write_text(
            new_tree.code, encoding="utf-8",)
    else:
        print(f"No changes needed in {target_file}")


if __name__ == "__main__":
    main()
