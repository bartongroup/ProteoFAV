#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
ProteFAV: protein feature aggregation and variants
--------------------------------------------------

Exploring the power of Pandas to work with protein structures,
sequences and genetic variants.

:copyright: (c) 2015-2017.
:license: TBD, see LICENSE for more details.
"""


from setuptools import setup
from setuptools import find_packages

from proteofav import __version__, __authors__


setup(
    # Basic package information.
    name='proteofav',
    version=__version__,
    packages=find_packages(),

    # # Packaging options.
    package_data={'': ['*.ipynb', '*.rst']},

    # Package dependencies.
    install_requires=['pandas', 'requests', 'scipy', 'numpy', 'lxml'],

    # Tests.
    test_suite='tests',

    # Test dependencies.
    test_requires=['mock;python_version<"3.4"'],

    # Metadata for PyPI.
    author=__authors__,
    author_email='tbrittoborges@dundee.ac.uk',
    license='TBD',
    url='https://github.com/bartongroup/ProteoFAV/tree/master',
    keywords='bioinformatics structural-biology data-analysis python pandas',
    description='PROtein Feature Aggregation and Variants.',
    long_description=open('README.rst').read(),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: TBD',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Programming Language :: Python :: Implementation :: PyPy3',
        'Topic :: Internet',
        'Topic :: Scientific/Engineering :: Bio-informatics',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]
)