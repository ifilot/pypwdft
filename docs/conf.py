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
import sys
from pathlib import Path

# Import the package from this checkout so autodoc always reflects the source
# being documented, irrespective of Sphinx's working directory or whether an
# older PyPWDFT release is installed in the build environment.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))
import sphinx_rtd_theme

# -- Project information -----------------------------------------------------

project = 'PyPWDFT'
copyright = '2024, Ivo Filot'
author = 'Ivo Filot'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.mathjax',
    'sphinx.ext.autosectionlabel',
    'sphinx_rtd_theme',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon'
]

autodoc_default_options = {
    'members': True,
    'show-inheritance': True,
}
autodoc_typehints = 'description'

# These runtime dependencies are not needed to inspect the public API, and
# documentation builders commonly run without the compiled scientific stack.
# Mocking them only during autodoc import keeps the reference sourced from the
# real PyPWDFT modules without requiring a working numerical backend.
autodoc_mock_imports = [
    'mendeleev',
    'pyfftw',
]

suppress_warnings = ['autosectionlabel.*']

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
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
master_doc = 'index'
html_static_path = ['_static']
# html_theme_options = {
#     'display_version': True,
#     'analytics_id': 'G-H71EPP6GVB'
# }
html_logo = "_static/img/pypwdft_128px.png"
html_favicon = "_static/img/favicon.ico"
html_css_files = [
    "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.1.1/css/all.min.css"
]

# other options
html_show_sourcelink = False

def setup(app):
   app.add_css_file('css/custom.css')
