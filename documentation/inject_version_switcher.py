"""Inject the FELiCS version switcher widget into every HTML page of a built
documentation version.

Runs once per version, right after that version's "sphinx-build" call (see the
"create documentation" job in .gitlab-ci.yml), so it also covers already-tagged
releases whose own historical conf.py predates this feature and cannot be
edited retroactively. The widget itself (felics_version_switcher.js/.css) is
published once, shared by every version, under public/_static/.
"""

import argparse
import os


HEAD_INJECTION_TEMPLATE = (
    '<link rel="stylesheet" href="{pages_url}/_static/felics_version_switcher.css">\n'
    '<script src="{pages_url}/_static/felics_version_switcher.js" defer></script>\n'
    '</head>'
)


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--pages-url',
        default='',
        help="Base URL the documentation is published under (CI_PAGES_URL); empty for a root-relative local build.",
    )
    parser.add_argument(
        'build_dir',
        help="Directory a single version was built into (e.g. public/v3.1.1).",
    )
    return parser.parse_args()


def inject_into_html_files(build_dir, pages_url):
    injection = HEAD_INJECTION_TEMPLATE.format(pages_url=pages_url)
    for root, _, file_names in os.walk(build_dir):
        for file_name in file_names:
            if not file_name.endswith('.html'):
                continue
            file_path = os.path.join(root, file_name)
            with open(file_path, 'r', encoding='utf-8') as html_file:
                content = html_file.read()
            if '</head>' not in content:
                continue
            with open(file_path, 'w', encoding='utf-8') as html_file:
                html_file.write(content.replace('</head>', injection, 1))


def main():
    arguments = parse_arguments()
    inject_into_html_files(arguments.build_dir, arguments.pages_url)


if __name__ == '__main__':
    main()
