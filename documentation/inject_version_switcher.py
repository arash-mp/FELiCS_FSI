"""Inject the FELiCS version switcher widget into every HTML page of a built
documentation version.

Runs once per version, right after that version's "sphinx-build" call (see the
"create documentation" job in .gitlab-ci.yml), so it also covers already-tagged
releases whose own historical conf.py predates this feature and cannot be
edited retroactively. The widget itself (felics_version_switcher.js/.css) is
published once, shared by every version, under public/_static/.

The injected <link>/<script> tags use a path relative to each individual page
(computed from its depth below build_dir, the same way Sphinx resolves its own
static assets) rather than an absolute CI_PAGES_URL, so the same build works
unmodified at the Pages site root, under a Pages path_prefix preview, or
opened straight from a downloaded artifact.
"""

import argparse
import os


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'build_dir',
        help="Directory a single version was built into (e.g. public/v3.1.1).",
    )
    return parser.parse_args()


def relative_static_prefix(build_dir, file_path):
    # felics_version_switcher.{js,css} are shared by every version and live
    # one level above build_dir (public/_static, not public/<version>/_static),
    # hence the "+ 1": even a page directly in build_dir has to step out of
    # its own version directory to reach it.
    depth = os.path.relpath(file_path, build_dir).count(os.sep)
    return '../' * (depth + 1)


def inject_into_html_files(build_dir):
    for root, _, file_names in os.walk(build_dir):
        for file_name in file_names:
            if not file_name.endswith('.html'):
                continue
            file_path = os.path.join(root, file_name)
            with open(file_path, 'r', encoding='utf-8') as html_file:
                content = html_file.read()
            if '</head>' not in content:
                continue
            static_prefix = relative_static_prefix(build_dir, file_path)
            injection = (
                f'<link rel="stylesheet" href="{static_prefix}_static/felics_version_switcher.css">\n'
                f'<script src="{static_prefix}_static/felics_version_switcher.js" defer></script>\n'
                '</head>'
            )
            with open(file_path, 'w', encoding='utf-8') as html_file:
                html_file.write(content.replace('</head>', injection, 1))


def main():
    arguments = parse_arguments()
    inject_into_html_files(arguments.build_dir)


if __name__ == '__main__':
    main()
