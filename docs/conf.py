import os
import sys

# Add project root to sys.path so autodoc can find miniBSE
sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------------
project = "miniBSE"
copyright = "2026, miniBSE Developers"
author = "Ivan Infante et al."
release = "1.0.0"
version = "1.0"

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx_copybutton",
    "myst_parser",
]

# Source suffix
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# Master document
master_doc = "index"

# Mock C++ extensions and optional heavy dependencies for doc building
autodoc_mock_imports = [
    "libint_cpp",
    "torch",
]

# Autodoc settings
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autoclass_content = "both"

# Napoleon settings (Google and NumPy docstrings)
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_use_param = True
napoleon_use_rtype = True

# Copybutton configuration
copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: "
copybutton_prompt_is_regexp = True

# Templates & static paths
templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

html_theme_options = {
    "collapse_navigation": False,
    "sticky_navigation": True,
    "navigation_depth": 4,
    "includehidden": True,
    "titles_only": False,
    "logo_only": False,
}

html_context = {
    "display_github": True,
    "github_user": "nlesc-nano",
    "github_repo": "miniBSE",
    "github_version": "main",
    "conf_py_path": "/docs/",
}
