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
import os
import sys
from sphinx_gallery.sorting import FileNameSortKey
import sphinx_gallery.backreferences

sphinx_gallery.backreferences.THUMBNAIL_TEMPLATE = (
    sphinx_gallery.backreferences.THUMBNAIL_TEMPLATE.replace("snippet", "title")
)

sys.path.insert(0, os.path.abspath(".."))


# -- Project information -----------------------------------------------------

project = "SimuPy Flight"
copyright = "2022, Benjamin W. L. Margolis"
author = "Benjamin W. L. Margolis"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.autodoc",
    "sphinx_gallery.gen_gallery",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

sphinx_gallery_conf = {
    "examples_dirs": os.path.join(os.path.dirname(__file__), "..", "nesc_test_cases"),
    "gallery_dirs": "nesc_test_cases",
    "filename_pattern": r"/nesc_case",
    "ignore_pattern": r"(F16_)|(nesc_testcase_helper)|(process_NESC)|(run_nesc).*\.py",
    "within_subsection_order": FileNameSortKey,
    "reset_argv": lambda gallery_conf, script_vars: ["--no-test"],
    "reference_url": {
        "simupy_flight": None,
        # "simupy": "https://simupy.readthedocs.io/en/latest/",
    },
}


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

html_css_files = ["custom.css"]
