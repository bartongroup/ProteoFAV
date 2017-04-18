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
to genetic variants. The tool relies heavily in `Pandas`_ library to quickly load protein data i
into DataFrames for data analysis in Python or exported and analysed elswhere. It excels
on integration of the protein features, structural data and genetic variation.

Installation and developing
~~~~~~~~~~~~~~~~~~~~~~~~~~~

With conda:
    conda-env create -n proteofav -f path/to/ProteoFAV/requirements.txt
    source activate proteofav
    cd path/to/ProteoFAV
    pip install .

for developing and testing. Test dependencies should be resolved with:
    python setup.py develop --user

Tests are run with:
    cd path/to/ProteoFAV/tests
    python -m unittest discover


Configuration
~~~~~~~~~~~~~

After installing run:
    proteofav-setup

To set-up the download directories for mmCIF (db_mmcif), SIFTS (db_sifts), DSSP (db_dssp),
Ensembl Germline (db_germline_variants) and Ensembl Somatic (db_somatic_variants) in the config.txt.


Usage
~~~~~

Working on documenting the package `docs`_...

Dependencies
~~~~~~~~~~~~

The framework was developed to support Python 2.7+ and Python 3.4+. Check
`requirements`_ for specific requirements.

Package architecture
~~~~~~~~~~~~~~~~~~~~

TBA

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
- Adding new file parsers and extending the functionality for user provider
- Extending parsers to support to edge cases

.. _requirements: https://github.com/bartongroup/ProteoFAV/blob/master/requirements.txt
.. _license: https://github.com/bartongroup/ProteoFAV/blob/master/LICENSE.txt
.. _issue tracker: https://github.com/bartongroup/ProteoFAV/issues
.. _docs: https://github.com/bartongroup/ProteoFAV/blob/master/docs/index.rst
.. _Pandas: http://pandas.pydata.org/