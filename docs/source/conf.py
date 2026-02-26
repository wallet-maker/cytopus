# Configuration file for the Sphinx documentation builder.

import os
import sys
sys.path.insert(0, os.path.abspath('..'))

project = 'Cytopus'
copyright = '2024'
author = 'Thomas Walle'
release = '2.0'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
]

templates_path = ['_templates']
exclude_patterns = []

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False,
}

html_static_path = ['_static']

def setup(app):
    app.add_css_file('custom.css')

# Intersphinx mapping for references
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
    'networkx': ('https://networkx.org/documentation/stable/', None),
}

# Path to a logo image from the repository `img/` folder (displayed in the RTD theme header)
# Relative path is from this file (docs/source)
html_logo = '../../img/cytopus_v1.1_stable_graph.png'

