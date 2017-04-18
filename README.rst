ProteoFAV
=========

**Protein feature aggregation and variants**


.. image:: https://img.shields.io/pypi/v/proteofav.svg
        :target: https://pypi.python.org/pypi/proteofav

.. image:: https://img.shields.io/travis/tbrittoborges/proteofav.svg
        :target: https://travis-ci.org/tbrittoborges/proteofav

.. image:: https://readthedocs.org/projects/proteofav/badge/?version=latest
        :target: https://proteofav.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://pyup.io/repos/github/tbrittoborges/proteofav/shield.svg
     :target: https://pyup.io/repos/github/tbrittoborges/proteofav/
     :alt: Updates


ProteFAV is a Python framework to fetch, process and integrate protein structure and features
to genetic variants. The tool relies heavily in `Pandas`_ library.

The tool aims to quickly load protein data into DataFrames for data analysis in Python.
Alternatively the data can be exported and analysed elsewhere. The tool excels in data
integration of the protein features, structural data and genetic variation.

Installation
~~~~~~~~~~~~

With conda:
    conda-env create -n proteofav -f path/to/ProteoFAV requirements.txt
    source activate proteofav
    cd path/to/ProteoFAV
    pip install .

Configuration
~~~~~~~~~~~~~

Set db_mmcif, db_sifts, db_dssp, db_germline_variants, db_somatic_variants with proteofav-setup,
 or simply modify config.txt. Set all to the download directory, such as
 /path/to/Downloads/proteofav or, for bigger projects, organise a project structure.


Testing
~~~~~~~

Test dependencies need to be installed prior running the tests. Currently ProteoFAV uses Unittest:

    cd path/to/Proteofav/tests
    python -m unittest discover

Usage
~~~~~

Working on documenting the package `docs`_...

Dependencies
~~~~~~~~~~~~

The framework was developed to support Python 2.7+ and Python 3.4+. Check
`requirements`_ for specific requirements.

Package architecture
~~~~~~~~~~~~~~~~~~~+

The framework's test and documentations are in /test and docs, respectevely.

- config: User configuration
- uniprot:
- variants:

Contributing and Bug tracking
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Feel free to fork, clone, share and distribute. If you find any bugs or
issues please log them in the `issue tracker`_.

License
~~~~~~~

See `license`_.

Future directions
~~~~~~~~~~~~~~~~~

We are currently working for several features and improvements:
- Using Python's object orientation to generalise the data integration, and so merge table
routine will get smarter
- Adding new file parsers and extending the functionality for user provider (instead public
avaible)
- Improving the parsers and extending the support to edge cases



.. _requirements: https://github.com/bartongroup/ProteoFAV/blob/master/requirements.txt
.. _license: https://github.com/bartongroup/ProteoFAV/blob/master/LICENSE.txt
.. _issue tracker: https://github.com/bartongroup/ProteoFAV/issues
.. _docs: https://github.com/bartongroup/ProteoFAV/blob/master/docs/index.rst
.. _Pandas: http://pandas.pydata.org/