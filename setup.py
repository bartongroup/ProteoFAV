# -*- coding: utf-8 -*-

"""
ProteoFAV: Protein Feature Aggregation and Variants
---------------------------------------------------

Open-source framework for simple and fast integration
of protein structure data with sequence annotations
and genetic variation

:copyright: (c) 2015-2017.
:license: TBD, see LICENSE for more details.
"""

from setuptools import setup, find_packages

from proteofav import __version__, __authors__


setup(
    # Basic package information.
    name='proteofav',
    version=__version__,
    packages=find_packages(),

    # Packaging options.
    # package_data={'': ['*.ipynb', '*.rst']},
    include_package_data=True,
    py_modules=['proteofav.main'],

    # Package dependencies.
    install_requires=['pandas>=0.17',
                      'requests>=2.12',
                      'lxml>=3.6',
                      'click>=6',
                      'scipy'
                      ],
    test_requires=['mock', 'python_version>"3.4"'], #

    # Tests.
    test_suite='tests',

    # Metadata for PyPI.
    author=__authors__,
    author_email='tbrittoborges@dundee.ac.uk',
    license='LICENSE.txt',
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