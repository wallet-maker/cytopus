import datetime
import os
import sys

# Add repo root so `import cytopus` resolves
sys.path.insert(0, os.path.abspath('../..'))

# -- Project information ---
project = 'cytopus'
copyright = f'{datetime.datetime.now().year}, Thomas Walle'
author = 'Thomas Walle'
try:
    from cytopus import __version__
    release = __version__
except ImportError:
    release = "2.0.0"

# -- General configuration ---
extensions = [
    'myst_parser',
    'sphinx.ext.autodoc',           # Auto-generate API docs
    'sphinx.ext.napoleon',          # Support for NumPy/Google docstring styles
    'sphinx.ext.viewcode',          # Link to source code
    'sphinx.ext.intersphinx',       # Link to other projects (e.g., numpy, scipy)
    'sphinx.ext.autosummary',       # Generate summary tables for modules/classes/functions
    'sphinx_design',                # Cards/grids for better layout
    'nbsphinx',                     # Jupyter notebook support
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']

# -- Autosummary configuration ---
autosummary_generate = True

# -- Autodoc configuration ---
autodoc_member_order = 'bysource'
autodoc_typehints = 'description'
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': False,
    'show-inheritance': True,
}
# -- Napoleon configuration (for docstring parsing) ---
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
# -- nbsphinx configuration ---
nbsphinx_execute = 'never'  # Don't execute notebooks on build (can set to 'always' later)
nbsphinx_kernel_name = 'python3'
# -- Intersphinx mapping ---
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
    'scanpy': ('https://scanpy.readthedocs.io/en/stable/', None),
    'anndata': ('https://anndata.readthedocs.io/en/stable/', None),
}

html_theme = 'furo'  # Modern, responsive theme
html_static_path = ['_static']
html_title = 'Cytopus Documentation'
html_css_files = ['output.css']

# -- Additional HTML context ---
html_context = {
    'github_user': 'wallet-maker',
    'github_repo': 'cytopus',
    'github_version': 'main',
    'doc_path': 'docs',
}