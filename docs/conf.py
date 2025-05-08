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

project = 'ctfr_benchmarks'
author = "Bernardo A. Boechat, Maurício do V. M. da Costa, Luiz W. P. Biscainho"
copyright = "2025, Bernardo A. Boechat, Maurício do V. M. da Costa, Luiz W. P. Biscainho"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.intersphinx"
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


# -- Generate code documentation -----------------------------------------------------

import os
from ctfr_bm import get_method_name

IMPLEMENTATIONS_BASE_PATH = "../src/ctfr/implementations"
CODE_FOLDER = "code"
CODE_TREE_FILE = "code_tree.rst"

def pyx_name_from_method_key(method_key):
    return f"{method_key}_cy.pyx"

methods_to_doc = [
    "baseline_swgm",
    "swgm",
    "baseline_fls",
    "fls",
    "baseline_lt",
    "lt",
    "baseline_sls",
    "sls_h",
    "sls_i"
]

code_tree_toctree_paths = []

for method in methods_to_doc:

    method_name = get_method_name(method)
    pyx_name = pyx_name_from_method_key(method)
    pyx_path = os.path.join(IMPLEMENTATIONS_BASE_PATH, pyx_name)
    code_output_path_no_extension = f"{os.path.join(CODE_FOLDER, method)}"
    code_tree_toctree_paths.append(code_output_path_no_extension)

    with open(code_output_path_no_extension + ".rst", "w") as f:
        title = f"{method_name}"
        f.write(title + "\n")
        f.write("="*len(title) + "\n\n")
        f.write(f"Contents of ``{pyx_name}``" + "\n\n")
        f.write(".. literalinclude:: " + "../" + pyx_path + "\n")

with open(CODE_TREE_FILE, "w") as f:
    f.write(".. toctree::" + "\n")
    f.write("   :caption: Source code (methods)" + "\n\n")
    for path in code_tree_toctree_paths:
        f.write(f"   {path}" + "\n")