"""Sphinx configuration for PyGAGE documentation."""
import os
import sys
from datetime import datetime

# Make the package importable for autodoc (repo root = two levels up).
sys.path.insert(0, os.path.abspath("../.."))

# -- Project information ------------------------------------------------------
project = "PyGAGE"
author = "Richard A. White III, Jose L. Figueroa III"
copyright = f"{datetime.now().year}, {author} (CC BY-NC 4.0)"

try:
    from pygage import __version__ as release
except Exception:  # pragma: no cover
    release = "1.2.0"
version = ".".join(release.split(".")[:2])

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",       # Google/NumPy-style docstrings
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "myst_parser",               # Markdown (CHANGELOG, etc.)
    "sphinx_copybutton",
]

templates_path = ["_templates"]
exclude_patterns = []
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
master_doc = "index"

# -- Autodoc / autosummary ---------------------------------------------------
autosummary_generate = True
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}
# Optional heavy deps that may be absent when building docs: mock so imports don't fail.
autodoc_mock_imports = ["anndata"]
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_use_rtype = False

# -- Intersphinx -------------------------------------------------------------
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "polars": ("https://docs.pola.rs/api/python/stable/", None),
}

# -- HTML output -------------------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "navigation_depth": 3,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "titles_only": False,
}
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_title = f"PyGAGE {version}"

# copybutton: strip prompts
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True

todo_include_todos = False
