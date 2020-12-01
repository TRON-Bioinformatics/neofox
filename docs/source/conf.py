# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = "NeoFox"
copyright = "2020, TRON – Translational Oncology at the University Medical Center of the Johannes Gutenberg University Mainz - Computational Medicine"
author = "Franziska Lang & Pablo Riesgo Ferreiro"

# The full version, including alpha/beta/rc tags
# TODO: it would be great to have a common versioning, but right now it would require importing neofox here which does
# TODO: not work in readthedocs environment
release = "0.4.0.dev2"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "m2r2",
    #'recommonmark'
    "nbsphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx.ext.autodoc",
]

source_suffix = [".rst", ".md"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'sphinx_rtd_theme'

# html_theme = 'insipid'
# html_theme_options = {
#    'body_centered': False,
#    'body_max_width': None,
#    'breadcrumbs': True,
#    'right_buttons': ['repo-buttom.html', 'tron-buttom.html']
# }
# html_context = {
#    'display_github': True,
#    'github_user': 'TRON-Bioinformatics',
#    'github_repo': 'neofox',
# }


html_theme = "pydata_sphinx_theme"
html_copy_source = False
# html_add_permalinks = '\N{SECTION SIGN}'


html_theme_options = {
    "github_url": "https://github.com/TRON-bioinformatics/neofox",
    "external_links": [
        {
            "name": "TRON",
            "url": "https://tron-mainz.de",
            "img": "_templates/tron-small.svg",
        }
    ],
}

html_logo = "_templates/tron.svg"
html_favicon = "_templates/tron_small.svg"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
