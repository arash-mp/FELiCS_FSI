# FELiCS developer guidelines for AI coding agents

FELiCS (*Finite Element Linearized Combustion Solver*) is a Python-based CFD tool for
linearized flow analysis (turbulence, heat/mass transport, chemical reactions, acoustics).

This file distills the developer guidelines from the
[project wiki](https://gitlab.com/felics-group/FELiCS/-/wikis/home) into instructions an AI
coding agent can follow directly. The wiki remains the canonical, human-maintained source —
if something here looks out of date or incomplete, check the wiki page linked in each section
and prefer it over this file.

## Required steps for every code change

Whenever you write new code or change existing code, before considering the work done:

1. Always run the unit tests: `python -m pytest` (see [Unit testing](#unit-testing) below).
2. Always run the full felics-tests suite (see [Test cases](#test-cases) below) — every case
   has to pass, not just the ones you think are related to your change. Read that section's
   gotchas first: a clean exit code from `runAll.sh` does **not** mean the tests passed.
3. Run the tutorials to confirm they still work: `bash tutorials/test_tutorials.sh` (requires an
   activated `felics` conda environment). This is separate from the `how_tos/` notebooks
   mentioned below.
4. Check whether your change requires the [documentation](#documentation) to be updated or new
   content to be added — new/changed public behavior, parameters, or governing equations
   usually need a corresponding doc update.

## Rules

- Never create a git commit without explicit human approval, even if all checks above pass.
- There is no GitLab CLI or API access in this environment. Never fabricate or silently skip a
  GitLab-side action (opening a merge request, assigning a reviewer) — hand it back to a human
  instead (see [Before finishing work](#before-finishing-work--opening-a-merge-request)).

## Coding conventions

Full guide: [wiki: Coding conventions](https://gitlab.com/felics-group/FELiCS/-/wikis/Coding-conventions).
Base style is [PEP 8](https://peps.python.org/pep-0008/), plus these FELiCS-specific rules:

- English only for names, comments, and docstrings.
- Prefer long, descriptive names over abbreviations (`number_grid_points`, not `n`).
- 4-space indentation, no tabs. One space around operators (`=`, `+`, `==`, `and`, ...), no
  space inside brackets, no space before `,`/`:`, one space after. No trailing whitespace.
- Trailing commas in multi-line calls/data structures (keeps Git diffs clean).
- Multi-line function signatures/calls: opening bracket ends the first line, one argument per
  line with a hanging indent, closing bracket on its own line.
- Imports: all at the top of the file, one per line, explicit (`from module import X, Y`,
  never `import *`), absolute rather than relative. Group as stdlib / third-party / local,
  each group separated by a blank line and a comment naming the group.
- Naming: classes in `PascalCase`, functions/methods in `snake_case`.
- Visibility: no leading underscore = public; single leading underscore (`_x`) = protected;
  double leading underscore (`__x`) = private (name-mangled). Expose internal attributes that
  need external read/write access via a `@property`, not by making them public.
- Within a class body, order methods private → protected → public, alphabetically within each
  group.
- Comments: don't restate the code. Prefer clear code that needs fewer comments over commented
  bad code. Comment things a reader might mistake as redundant (workarounds, non-obvious
  assumptions, perf tricks). Prefix adaptation notes with `TODO:`. Leave a blank line before a
  comment.
- Docstrings: [numpydoc format](https://numpydoc.readthedocs.io/en/latest/format.html)
  ([examples](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html)).
  - Every class starts with a one-line summary, optionally followed by a longer description.
  - `__init__` docstrings are not rendered (it's treated as private), so document constructor
    parameters in the *class* docstring instead, under its own `Parameters` section.
  - List public methods under a `Methods` section in the class docstring — but do not
    duplicate that method list inside the class's own top docstring block if it's the same
    one already documented per-method.

## Test cases

FELiCS-tests: the [felics-tests](https://gitlab.com/felics-group/felics-tests) repo holds the
end-to-end validation cases for FELiCS.

- **VERY IMPORTANT:** To validate new or changed code, run the **full** felics-tests suite
  locally against your changes before considering the work done — every case has to pass, not
  just the ones you think are related. This is required in addition to the unit tests below,
  not a replacement for them: unit tests catch regressions at the function/class level,
  felics-tests validates end-to-end behavior. Do not skip either.

  How to run felics-tests locally:
  1. Clone the repo separately from FELiCS: `git clone https://gitlab.com/felics-group/felics-tests`.
  2. Check out its `development` branch and merge/rebase in the latest `development` if needed.
  3. Make sure FELiCS itself (this repo, with your changes) is installed correctly into your
     Python environment.
  4. Activate the `felics` conda environment.
  5. From the top level of the felics-tests repo, run `source runAll.sh` (see gotchas below —
     do not run it as `./runAll.sh`).
  6. Wait for every case to report a result, then check that all say `successful` (see below).
     Expect this to take a few minutes — don't interrupt it early or assume a hang.

  Three gotchas that trip up non-interactive/agent shells specifically:
  - **The `FELiCS` command is a shell function, not a program on `PATH`.** It's defined in
    `~/.bashrc` (wraps `python .../src/main.py`). A plain non-interactive shell (e.g. `bash -c
    "..."`, or most AI-agent shell tools) does not source `~/.bashrc`, so `FELiCS` will fail
    with `command not found` unless you explicitly `source ~/.bashrc` (or otherwise define the
    function) first in that same shell before calling `runAll.sh`.
  - **`runAll.sh` is not marked executable and must be sourced, not run as a script.** It relies
    on `FELiCS` being a function in the *current* shell, so it only works via `source runAll.sh`
    (or `. runAll.sh`) in a shell that already has the function defined — running it as
    `./runAll.sh` or `bash runAll.sh` spawns a fresh subshell that never sees the function, again
    producing `command not found`.
  - **`runAll.sh` backgrounds every case with `&` and never `wait`s** — it returns almost
    immediately while the actual FELiCS runs and comparisons continue in the background. A
    non-zero/zero exit code or the script "finishing" tells you nothing. Instead, wait until all
    child processes are done (e.g. poll for `python compareResults.py` processes under your
    shell), then check the captured output for one `Validation case N/12 '<name>' successful!
    Maximum residuum is ...` line per case (12 total as of this writing — the count is in each
    printed line). Any line saying `failed!` instead of `successful!` is a real failure, not
    just non-zero shell status.

## Unit testing

Full guide: [wiki: Unit test](https://gitlab.com/felics-group/FELiCS/-/wikis/unit-test).

- Any new functionality or method needs a corresponding unit test.
- When modifying existing functionality, run the related unit tests and update them if they
  fail (after confirming the failure is expected, not a regression).
- Two accepted test concepts: **validation** (compare against an analytical reference
  implemented in the test) and **snapshot** (compare against a recorded prior output).
- Run the full suite from the repo root: `python -m pytest` (see [tests/README.md](tests/README.md)).
- These tests also run automatically in CI on pushes to `development`.


## Before finishing work / opening a merge request

Full checklist: [wiki: Checklist for merging into development](https://gitlab.com/felics-group/FELiCS/-/wikis/Workflow-merging-into-development).

**Assume there is no GitLab CLI or API access available.** Steps 1–5 below (git operations,
running checks locally) are things an AI coding agent can do directly. Steps 6–7 involve
GitLab-side actions (opening a merge request, assigning a reviewer) that an agent without
GitLab credentials cannot perform — do not attempt to fake or skip these. Instead, finish
steps 1–5, then hand the branch back to the human with a short note on what was done and
which of steps 6–7 still need to be done manually.

1. Merge the latest `development` into your branch.
2. Run the style checker from the repo root: `python tools/code_refactoring/felics_check_code.py src/`
   and make sure the code conforms to the coding conventions above.
3. Confirm the [required steps for every code change](#required-steps-for-every-code-change)
   (unit tests, full felics-tests suite, tutorials) are all still green.
4. If you touched anything used for scripting, also run the notebooks in `how_tos/` (separate
   from the `tutorials/` check above) and confirm they still execute end to end.
5. If you touched `documentation/` or docstrings, build the docs locally (see below) and check
   the rendered output.
6. Merge requests target the `development` branch, not `main` (`main` only receives merges at
   release time — see [wiki: FELiCS Release Guide](https://gitlab.com/felics-group/FELiCS/-/wikis/Workflows/FELiCS-Release-Guide)).
   Assign yourself, and **always assign at least one other person as reviewer** — never merge
   without a second person reviewing (current reviewers per the wiki's
   [List of people](https://gitlab.com/felics-group/FELiCS/-/wikis/List-of-people): Sophie,
   Thomas, Simon — verify against that page, it changes). *(GitLab-side step — see note above.)*
7. If there's a companion test case in the [felics-tests](https://gitlab.com/felics-group/felics-tests)
   repo, make sure its branch also has the latest `development` merged in and meets that repo's
   test guidelines. *(GitLab-side step — see note above.)*

## Branch naming

Full guide: [wiki: Branches: naming conventions](https://gitlab.com/felics-group/FELiCS/-/wikis/conventions-branch-names).

- `<issue-number>-<title>` — created automatically by GitLab from an issue/merge request.
- `dev_...` — a feature or development-branch change, merged into `development` and deleted
  when done.
- `personal_<name>_...` — personal/experimental branches, never merged upstream.
- `script_...` — branches tracking user scripts for compatibility.

## Documentation

Full guide: [wiki: Documentation Creation Guide](https://gitlab.com/felics-group/FELiCS/-/wikis/Documentation-Creation-Guide-for-Developers).

- Docs and doc-string changes are committed to `development`, with the commit message prefixed
  `[doc]`.
- Guides live under `documentation/` in three sections: `GoverningEquations`, `Tutorials`,
  `How-To-Guides`. New pages follow the template at `documentation/markdown_template.md` and
  must be added to the relevant section's `index` file.
- Docs are built with Sphinx + myst-parser + sphinx-book-theme + sphinx-autoapi. Build locally
  with `sphinx-build -b html . _build` from `documentation/`, then open `_build/index.html`.
- Use MyST admonitions (`` ```{note} ``, `` ```{warning} ``, `` ```{tip} ``, etc.) for callouts —
  see [wiki: Write nice markdown pages](https://gitlab.com/felics-group/FELiCS/-/wikis/Write-nice-markdown-pages-for-the-docu)
  for the full list.

## Issues and reviews

Full guide: [wiki: Milestones, issuing and merge requests](https://gitlab.com/felics-group/FELiCS/-/wikis/Milestones,-issuing-and-merge-requests).

- Work should be tied to a GitLab issue with a descriptive name, the working branch mentioned
  in the description, self-assigned, and a time estimate.
- Merge requests are best created directly from the issue so the branch/issue link is
  automatic.
