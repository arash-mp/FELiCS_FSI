#!/usr/bin/env python3
"""
FELiCS style checker (CI-friendly)

Checks (no rewriting):
  1) Line length > 80 that cannot be fixed by your rule:
        - only split if an eligible comma exists (comma token inside (),[],{})
     => report only lines >80 with NO eligible comma.
  2) Class names PascalCase (excluding EXCL_RENAME)
  3) Function names snake_case (excluding EXCL_RENAME)
  4) Imports must be top-level (excluding exact EXCL_IMPORTS)

Exit codes:
  0 = clean
  2 = violations found
"""

from __future__ import annotations

import io
import pathlib
import re
import sys
import tokenize
from dataclasses import dataclass
from typing import Dict, List

import libcst as cst
from libcst import metadata


# ---- keep these in sync with your rewriter ----
EXCL_RENAME = {"t", "getKSP", "getSize", "getVecs"}

EXCL_IMPORTS = {
    "iPython",
    "from fenics import project",
    "from FELiCS.Fields.MeanFlowClass import MeanFlowVertexValues",
    "from FELiCS.Fields.MeanFlowClass import MeanFlowClass",
    "from FELiCS.SpaceDisc.FEMSpaces import FEMSpaces",
}


def _norm_import_line(s: str) -> str:
    return " ".join(s.strip().split()).lower()


EXCL_IMPORTS_NORM = {_norm_import_line(x) for x in EXCL_IMPORTS}


@dataclass(frozen=True)
class Violation:
    file: str
    line: int
    col: int
    rule: str
    message: str
    snippet: str


def _is_pascal(name: str) -> bool:
    return bool(name) and name[0].isupper() and "_" not in name


def _is_snake(name: str) -> bool:
    return bool(re.fullmatch(r"[a-z_][a-z0-9_]*", name))


def _eligible_commas_by_line(source: str) -> Dict[int, List[int]]:
    """
    line -> comma columns for commas that are real tokens AND inside brackets.
    """
    toks = list(tokenize.generate_tokens(io.StringIO(source).readline))
    bracket_depth = 0
    out: Dict[int, List[int]] = {}

    for ttype, tstr, (srow, scol), _end, _line in toks:
        if ttype != tokenize.OP:
            continue

        if tstr in "([{":
            bracket_depth += 1
        elif tstr in ")]}":
            bracket_depth = max(0, bracket_depth - 1)
        elif tstr == "," and bracket_depth > 0:
            out.setdefault(srow, []).append(scol)

    return out


class StyleReportVisitor(cst.CSTVisitor):
    METADATA_DEPENDENCIES = (metadata.PositionProvider,
                             metadata.ParentNodeProvider)

    def __init__(self, file_path: str, source_lines: List[str]) -> None:
        self.file_path = file_path
        self.source_lines = source_lines
        self.violations: List[Violation] = []

    def _add(self, node: cst.CSTNode, rule: str, message: str) -> None:
        pos = self.get_metadata(metadata.PositionProvider, node)
        line = pos.start.line
        col = pos.start.column
        snippet = self.source_lines[line -
                                    1] if 1 <= line <= len(self.source_lines) else ""
        self.violations.append(
            Violation(
                file=self.file_path,
                line=line,
                col=col,
                rule=rule,
                message=message,
                snippet=snippet,
            )
        )

    def visit_ClassDef(self, node: cst.ClassDef) -> None:
        name = node.name.value
        if name in EXCL_RENAME:
            return
        if not _is_pascal(name):
            self._add(node.name, "class-name",
                      f"Class '{name}' is not PascalCase")

    def visit_FunctionDef(self, node: cst.FunctionDef) -> None:
        name = node.name.value
        if name in EXCL_RENAME:
            return
        if not _is_snake(name):
            self._add(node.name, "function-name",
                      f"Function '{name}' is not snake_case")

    def visit_SimpleStatementLine(self, node: cst.SimpleStatementLine) -> None:
        has_import = any(isinstance(x, (cst.Import, cst.ImportFrom))
                         for x in node.body)
        if not has_import:
            return

        norm = _norm_import_line(cst.Module(body=(node,)).code)
        if norm in EXCL_IMPORTS_NORM:
            return

        parent = self.get_metadata(metadata.ParentNodeProvider, node)
        if not isinstance(parent, cst.Module):
            self._add(node, "import-placement",
                      "Import is not at top-level (nested import)")


def _group_report_by_rule(violations: List[Violation]) -> str:
    """
    Grouped report:
      RULE
        file:line:col message
          snippet
    """
    by_rule: Dict[str, List[Violation]] = {}
    for v in violations:
        by_rule.setdefault(v.rule, []).append(v)

    out: List[str] = []
    for rule in sorted(by_rule):
        vs = sorted(by_rule[rule], key=lambda x: (x.file, x.line, x.col))
        out.append(f"\n[{rule}] ({len(vs)})")
        out.append("-" * (len(rule) + 4 + len(str(len(vs))) + 3))

        current_file = None
        for v in vs:
            if v.file != current_file:
                current_file = v.file
                out.append(f"\n  {current_file}")

            out.append(f"    - {v.line}:{v.col}  {v.message}")
            if v.snippet.strip():
                out.append(f"      {v.snippet.rstrip()}")
    return "\n".join(out).lstrip("\n")


def check_tree(root: pathlib.Path, max_len: int = 80) -> List[Violation]:
    violations: List[Violation] = []

    for pyfile in root.rglob("*.py"):
        try:
            source = pyfile.read_text(encoding="utf-8")
        except Exception:
            continue

        source_lines = source.splitlines()
        eligible = _eligible_commas_by_line(source)

       # Rule: line length > 80 AND HAS an eligible comma => violation
        for i, line in enumerate(source_lines, start=1):
            if len(line) <= max_len:
                continue
            if eligible.get(i):  # only flag if it CAN be split by your rule
                violations.append(
                    Violation(
                        file=str(pyfile),
                        line=i,
                        col=max_len,
                        rule="line-length",
                        message=f"Line > {max_len} chars and has an eligible comma (should be split)",
                        snippet=line,
                    )
                )

        # CST-based checks
        try:
            module = cst.parse_module(source)
        except Exception:
            # parse errors: treat as a violation (optional). Here: report and continue.
            violations.append(
                Violation(
                    file=str(pyfile),
                    line=1,
                    col=0,
                    rule="parse-error",
                    message="Failed to parse file with LibCST",
                    snippet="",
                )
            )
            continue

        wrapper = metadata.MetadataWrapper(module)
        visitor = StyleReportVisitor(str(pyfile), source_lines)
        wrapper.visit(visitor)
        violations.extend(visitor.violations)

    return violations


def main() -> None:
    if len(sys.argv) != 2:
        print("Usage: felics_style_check.py <target_dir>")
        sys.exit(1)

    root = pathlib.Path(sys.argv[1]).resolve()
    violations = check_tree(root, max_len=80)

    if not violations:
        print("✅ FELiCS style check passed.")
        sys.exit(0)

    print("❌ FELiCS style check failed.")
    print(_group_report_by_rule(violations))
    print(f"\nSummary: {len(violations)} total issues.")
    sys.exit(2)


if __name__ == "__main__":
    main()
