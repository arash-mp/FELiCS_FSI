#!/usr/bin/env python3
"""
Rope-based project-wide refactoring:
- Reads felics_renames.json (written by felics_style_rewriter.py).
- Applies renames across the project so imports/usages match new names.
"""

# Standard libraries
import json
import pathlib
import sys

# Third party libraries
from rope.base.project import Project
from rope.refactor.rename import Rename


def apply_rename(
    project: Project,resource_path: str,old_name: str,new_name: str,):
    resource = project.get_resource(resource_path)
    renamer = Rename(
        project,resource,old_name,)
    changes = renamer.get_changes(new_name)
    project.do(changes)


def main(
    root: pathlib.Path,):
    project = Project(str(root))
    log_file = root / "felics_renames.json"
    if not log_file.exists():
        print("No renames log found. Skipping rope refactor.")
        return

    try:
        renames = json.loads(log_file.read_text(encoding="utf-8"))
    except Exception as e:
        print(f"Failed to read renames log: {e}")
        return

    # Renames is a list of [file_path, old_name, new_name]
    applied = 0
    for entry in renames:
        try:
            resource_path, old_name, new_name = entry
            apply_rename(
                project,resource_path,old_name,new_name,)
            print(f"Rope: {resource_path}: {old_name} -> {new_name}")
            applied += 1
        except Exception as e:
            print(f"Rope: skip {entry}: {e}")

    print(f"Rope applied {applied} renames.")


if __name__ == "__main__":
    target = pathlib.Path(sys.argv[1]).resolve()
    main(target)
