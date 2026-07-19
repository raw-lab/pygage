"""Sphinx configuration for PyGAGE documentation."""

import os
import sys
from datetime import datetime

# Make the package importable for autodoc (src/ layout).
sys.path.insert(0, os.path.abspath("../../src"))

# -- Project information ------------------------------------------------------
project = "PyGAGE"
author = "Richard A. White III, Jose L. Figueroa III"
copyright = f"{datetime.now().year}, {author} \u00b7 CC BY-NC 4.0"

try:
    from pygage import __version__ as release
except Exception:  # pragma: no cover
    release = "1.2.1"
version = ".".join(release.split(".")[:2])

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",       # Google/NumPy-style docstrings
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "myst_parser",               # Markdown (CHANGELOG, etc.)
    "sphinx_copybutton",
    "sphinx_design",             # grids, cards, badges (landing page)
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

# -- HTML output (Furo) ------------------------------------------------------
html_theme = "furo"
html_title = f"PyGAGE {version}"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_logo = "_static/logo.svg"
html_favicon = "_static/favicon.svg"
html_copy_source = False
html_show_sphinx = False
pygments_style = "friendly"
pygments_dark_style = "material"

_BRAND = "#6d5efc"
_BRAND_DARK = "#a99cff"

html_theme_options = {
    "sidebar_hide_name": False,
    "navigation_with_keys": True,
    "light_css_variables": {
        "color-brand-primary": _BRAND,
        "color-brand-content": _BRAND,
        "color-api-name": _BRAND,
        "color-api-pre-name": _BRAND,
        "font-stack": ("Inter, -apple-system, BlinkMacSystemFont, 'Segoe UI', "
                       "Roboto, Helvetica, Arial, sans-serif"),
        "font-stack--monospace": ("'JetBrains Mono', 'SFMono-Regular', Menlo, "
                                  "Consolas, monospace"),
    },
    "dark_css_variables": {
        "color-brand-primary": _BRAND_DARK,
        "color-brand-content": _BRAND_DARK,
        "color-api-name": _BRAND_DARK,
        "color-api-pre-name": _BRAND_DARK,
    },
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/raw-lab/pygage",
            "html": (
                '<svg stroke="currentColor" fill="currentColor" viewBox="0 0 16 16">'
                '<path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38'
                "0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13"
                "-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66"
                ".07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15"
                "-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0"
                "1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82"
                "1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01"
                '1.93-.01 2.2 0 .21.15.46.55.38A8.01 8.01 0 0 0 16 8c0-4.42-3.58-8-8-8z">'
                "</path></svg>"
            ),
            "class": "",
        },
    ],
}

# copybutton: strip prompts so examples copy cleanly
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True
