# -*- coding: utf-8 -*-

"""
ProteoFAV: Protein Feature Aggregation and Variants
---------------------------------------------------

Open-source framework for simple and fast integration
of protein structure data with sequence annotations
and genetic variation

:copyright: (c) 2015-2017.
:license: MIT, see LICENSE for more details.
"""

import os
from setuptools import setup, find_packages

from proteofav import __version__, __authors__


def gather_dependencies():
    with open('requirements.txt', 'r') as f_in:
        return [l for l in f_in.read().rsplit(os.linesep)
                if l and not l.startswith("#")]
DEPENDENCIES = gather_dependencies()


setup(
    # Basic package information.
    name='proteofav',
    version=__version__,
    packages=find_packages(),

    # Packaging options.
    # package_data={'': ['*.ipynb', '*.rst']},
    package_data={'proteofav': ['*.ini']},
    include_package_data=True,
    py_modules=['proteofav.main'],

    # Package dependencies.
    install_requires=DEPENDENCIES,
    test_requires=['mock', 'python_version>"3.5"'],

    # Tests.
    test_suite='tests',

    # Metadata for PyPI.
    author=__authors__,
    author_email='tbrittoborges@dundee.ac.uk',
    license='LICENSE.md',
    entry_points={
        "console_scripts": ["proteofav-setup=proteofav.main:setup",
                            "proteofav=proteofav.main:main"]
        },
    url='https://github.com/bartongroup/ProteoFAV/tree/master',
    download_url="https://github.com/bartongroup/ProteoFAV/archive/master.zip",
    keywords='bioinformatics structural-biology data-analysis python pandas',
    description='PROtein Feature Aggregation and Variants.',
    long_description=open('README.rst').read(),

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Topic :: Internet',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'License :: OSI Approved :: MIT License',
    ]
)