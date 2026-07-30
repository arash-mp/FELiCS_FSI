"""Generate the version switcher JSON and the root redirect page for the
multi-version FELiCS documentation deployed to GitLab Pages.

Run once, after every version has been built into its own "public/<version>"
subdirectory (see the "create documentation" job in .gitlab-ci.yml). URLs are
deliberately relative to this directory rather than an absolute CI_PAGES_URL,
so the same output works unmodified at the Pages site root, under a Pages
path_prefix preview, or opened straight from a downloaded artifact.
"""

import argparse
import json
import os


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--preferred',
        required=True,
        help="Version identifier the site root should redirect to and mark as preferred.",
    )
    parser.add_argument(
        '--output-dir',
        required=True,
        help="Directory to write switcher.json and index.html into.",
    )
    parser.add_argument(
        'versions',
        nargs='+',
        help="Version identifiers that were built (e.g. development v3.1.1 v3.1.0).",
    )
    return parser.parse_args()


def build_switcher_entries(versions, preferred):
    def sort_key(version):
        if version == preferred:
            return (0, version)
        if version == 'development':
            return (2, version)
        return (1, version)

    entries = []
    for version in sorted(versions, key=sort_key):
        name = 'development (latest)' if version == 'development' else version
        entry = {
            'name': name,
            'version': version,
            'url': f'{version}/',
        }
        if version == preferred:
            entry['preferred'] = True
        entries.append(entry)
    return entries


def write_redirect_index(output_dir, preferred):
    redirect_url = f'{preferred}/index.html'
    content = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta http-equiv="refresh" content="0; url={redirect_url}">
  <link rel="canonical" href="{redirect_url}">
  <title>FELiCS documentation</title>
</head>
<body>
  <p>Redirecting to the <a href="{redirect_url}">{preferred}</a> documentation.</p>
</body>
</html>
"""
    with open(os.path.join(output_dir, 'index.html'), 'w') as index_file:
        index_file.write(content)


def write_switcher_json(output_dir, entries):
    with open(os.path.join(output_dir, 'switcher.json'), 'w') as switcher_file:
        json.dump(entries, switcher_file, indent=2)


def main():
    arguments = parse_arguments()
    entries = build_switcher_entries(arguments.versions, arguments.preferred)
    write_switcher_json(arguments.output_dir, entries)
    write_redirect_index(arguments.output_dir, arguments.preferred)


if __name__ == '__main__':
    main()
