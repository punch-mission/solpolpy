#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Note: To use the 'upload' functionality of this file, you must:
#   $ pipenv install twine --dev

import io
import os
import sys
#from numpy import get_include as np_get_include
#import numpy
from shutil import rmtree

from setuptools import find_packages, setup, Command
from setuptools.extension import Extension


# Package meta-data.
NAME = 'solpolpy'
DESCRIPTION = 'General-purpose solar depolarization/polarization suite'
URL = ''
EMAIL = 'mwest@boulder.swri.edu'
AUTHOR = 'Matthew J West, Marcus Hughes, Bryce Walbridge'
REQUIRES_PYTHON = '>=3.9.0'
VERSION = '0.0.1'

# What packages are required for this module to be executed?
REQUIRED = [
    'pytest',     # testing
    'numpy',      # basic numerics
    'astropy',    # FITS handling; WCS object
]

# What packages are optional?
EXTRAS = {
        'sunpy':    'sunpy',      # Required to be able to export SunPy objects
        'plotting': 'matplotlib',
        'docs':     'jupyter',    # to read the ipynb doc file
}

# The rest you shouldn't have to touch too much :)
# ------------------------------------------------
# Except, perhaps the License and Trove Classifiers!
# If you do change the License, remember to change the Trove Classifier for that!

here = os.path.abspath(os.path.dirname(__file__))

# Import the README and use it as the long-description.
# Note: this will only work if 'README.md' is present in your MANIFEST.in file!
try:
    with io.open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION

# Load the package's __version__.py module as a dictionary.
about = {}
if not VERSION:
    project_slug = NAME.lower().replace("-", "_").replace(" ", "_")
    with open(os.path.join(here, project_slug, '__version__.py')) as f:
        exec(f.read(), about)
else:
    about['__version__'] = VERSION

# Set up extensions - trivial with just helpers, but 
# useful later if we have to link in C libraries etc.
extensions = []

# Where the magic happens:
setup(
    name=NAME,
    version=about['__version__'],
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=find_packages(where='solpolpy'),
    # If your package is a single module, use this instead of 'packages':
    # py_modules=['mypackage'],

    # entry_points={
    #     'console_scripts': ['mycli=mymodule:cli'],
    # },
    requires=REQUIRED,
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    include_package_data=True,
    license='LGPL-v3',
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
    ],
    zip_safe = False,
)