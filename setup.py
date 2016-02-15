#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
ProteFAV: protein feature aggregation and variants
--------------------------------------------------

Exploring the power of Pandas to work with protein structures,
sequences and genetic variants.

:copyright: (c) 2015-2016.
:license: TBD, see LICENSE for more details.
"""


from setuptools import setup
from setuptools import find_packages

from proteofav import __version__


setup(
    # Basic package information.
    name='proteofav',
    version=__version__,
    packages=find_packages(),

    # Packaging options.
    include_package_data=True,

    # Package dependencies.
    install_requires=['pandas', 'requests', 'scipy', 'numpy', 'lxml', 'biopython'],

    # Tests.
    test_suite='tests',

    # Test dependencies.
    test_require=['mock', 'responses'],

    # Metadata for PyPI.
    # author=__author__,
    # author_email=__email__,
    license='MIT',
    url='http://github.com/biomadeira/proteofav/tree/master',
    keywords='python pandas pdb structures variants',
    description='Protein feature aggregation and variants.',
    long_description=open('README.rst').read(),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Programming Language :: Python :: Implementation :: PyPy3',
        'Topic :: Internet',
        'Topic :: Scientific/Engineering :: Bio-informatics',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]
)