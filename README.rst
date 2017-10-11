ProteoFAV
=========

**Protein Feature Aggregation and Variants**


.. image:: https://img.shields.io/pypi/v/proteofav.svg
        :target: https://pypi.python.org/pypi/proteofav

.. image:: https://img.shields.io/travis/bartongroup/proteofav.svg
        :target: https://travis-ci.org/bartongroup/proteofav

.. image:: https://readthedocs.org/projects/proteofav/badge/?version=latest
        :target: https://proteofav.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://pyup.io/repos/github/bartongroup/proteofav/shield.svg
     :target: https://pyup.io/repos/github/bartongroup/proteofav/
     :alt: Updates

ProteoFAV is a Python module that address the challenge of cross-mapping protein structures and protein sequences,
allowing for protein structures to be annotated with sequence features. It implements methods for working with
protein structures (via mmCIF, PDB, PDB Validation, DSSP and SIFTS files), sequence Features (via UniProt GFF annotations) and
genetic variants (via UniProt/EBI Proteins API and Ensembl REST API). Cross-mapping of structure and sequence is
performed with the aid of SIFTS.

ProteFAV relies heavily in the `Pandas`_ library to quickly load data into DataFrames for fast
data exploration and analysis. Structure and sequence
data are parsed/fetched onto Pandas DataFrames that are then merged-together (collapsed) onto a
single DataFrame.


Dependencies
~~~~~~~~~~~~

The framework was developed to support Python 3.5+ and Pandas 0.20+.

Check `requirements`_ for specific requirements.


Installation and developing
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Getting ProteoFAV:

.. code-block:: bash

    $ wget https://github.com/bartongroup/ProteoFAV/archive/master.zip -O ProteoFAV.zip
    $ unzip ProteoFAV.zip

    # alternatively
    $ git clone https://github.com/bartongroup/ProteoFAV.git


Installing With conda:

.. code-block:: bash

    $ conda-env create -n proteofav -f path/to/ProteoFAV/requirements.txt
    $ source activate proteofav
    $ cd path/to/ProteoFAV
    $ pip install .

Installing with Virtualenv:

.. code-block:: bash

    $ virtualenv --python `which python` env
    $ source env/bin/activate
    $ pip install -r requirements.txt


Test dependencies should be resolved with:

.. code-block:: bash

    $ python path/to/ProteoFAV/setup.py develop --user


Run the Tests with:

.. code-block:: bash

    $ python path/to/ProteoFAV/setup.py test

or:

.. code-block:: bash

    $ cd path/to/ProteoFAV/tests
    $ python -m unittest discover


Configuration
~~~~~~~~~~~~~

After installing run:

.. code-block:: bash

    $ proteofav-setup

To set-up the download directories for mmCIF (`db_mmcif`), SIFTS (`db_sifts`), DSSP (`db_dssp`),
PDB Validation (db_validation) and Annotations (db_annotation) in the
`config.ini`, otherwise ProteoFAV will download files to temporary directories.

Usage
~~~~~

TODO


Contributing and Bug tracking
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Feel free to fork, clone, share and distribute. If you find any bugs or
issues please log them in the `issue tracker`_.

Before you submit your *Pull-requests* read the `Contributing Guide`_.


Changelog
~~~~~~~~~

See the `Changelog`_

Licensing
~~~~~~~~~

See `LICENSE`_.

Credits
~~~~~~~

See the `Credits`_

.. _requirements: https://github.com/bartongroup/ProteoFAV/blob/master/requirements.txt
.. _LICENSE: https://github.com/bartongroup/ProteoFAV/blob/master/LICENSE.md
.. _issue tracker: https://github.com/bartongroup/ProteoFAV/issues
.. _docs: https://github.com/bartongroup/ProteoFAV/blob/master/docs/index.rst
.. _Pandas: http://pandas.pydata.org/
.. _Contributing Guide: https://github.com/bartongroup/ProteoFAV/wiki/Contributing-Guide
.. _Changelog: https://github.com/bartongroup/ProteoFAV/blob/master/CHANGELOG.rst
.. _Credits: https://github.com/bartongroup/ProteoFAV/blob/master/AUTHORS.rst
