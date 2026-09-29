"""Fail if a Sphinx HTML build shows raw RST math roles to readers.

Usage: python audit/check_rendered_math.py /path/to/sphinx/html
Sphinx's warning-free build alone does not catch roles inside strong markup.
"""
from pathlib import Path
import re
import sys

root = Path(sys.argv[1])
bad = []
for page in root.rglob("*.html"):
    if "_modules" in page.parts:  # autodoc exposes literal source code
        continue
    html = page.read_text(errors="replace")
    html = re.sub(r"<pre\b[^>]*>.*?</pre>", "", html, flags=re.I | re.S)
    html = re.sub(r"<code\b[^>]*>.*?</code>", "", html, flags=re.I | re.S)
    count = html.count(":math:")
    if count:
        bad.append((page.relative_to(root), count))
if bad:
    for page, count in bad:
        print(f"{page}: {count} raw :math: role(s)")
    raise SystemExit(1)
print("No raw :math: roles in reader-facing HTML pages.")
