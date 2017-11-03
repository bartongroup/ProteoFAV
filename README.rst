ProteoFAV
=========

**Protein Feature Aggregation and Variants**

.. |Pypi| |Build Status| |Documentation| |Coverage Status| |Health| |Pyup| |License|

|Build Status| |Documentation| |Python: versions| |License|

.. |Pypi| image:: https://img.shields.io/pypi/v/proteofav.svg
  :target: https://pypi.python.org/pypi/proteofav
.. |Build Status| image:: https://img.shields.io/travis/bartongroup/proteofav.svg
  :target: https://travis-ci.org/bartongroup/proteofav
.. |Documentation| image:: https://readthedocs.org/projects/proteofav/badge/?version=latest
  :target: https://proteofav.readthedocs.io/en/latest/?badge=latest
  :alt: Documentation Status
.. |Coverage Status| image:: https://coveralls.io/repos/github/bartongroup/proteofav/badge.svg?branch=master
  :target: https://coveralls.io/github/bartongroup/proteofav?branch=master
.. |Health| image:: https://landscape.io/github/bartongroup/proteofav/master/landscape.svg?style=flat
  :target: https://landscape.io/github/bartongroup/proteofav/master
.. |Pyup| image:: https://pyup.io/repos/github/bartongroup/proteofav/shield.svg
   :target: https://pyup.io/repos/github/bartongroup/proteofav/
   :alt: Updates
.. |License| image:: http://img.shields.io/badge/license-GPLv3-brightgreen.svg?style=flat
  :target: https://github.com/bartongroup/proteofav//blob/master/LICENSE.md
.. |Python: versions| image:: https://img.shields.io/badge/python-3.5,_3.6-blue.svg?style=flat
   :target: http://travis-ci.org/bartongroup/proteofav/

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
    # or
    $ cd path/to/ProteoFAV/tests
    $ python -m unittest discover


Configuration
~~~~~~~~~~~~~

After installing run:

.. code-block:: bash

    $ proteofav-setup

To set-up the download directories for mmCIF (``db_mmcif``), SIFTS (``db_sifts``), DSSP (``db_dssp``),
PDB Validation (``db_validation``) and Annotations (``db_annotation``) in the
``config.ini``, otherwise ProteoFAV will download files to temporary directories.

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

Credits
~~~~~~~

See the `Credits`_

Licensing
~~~~~~~~~

GNU General Public License v3 (GPLv3). See `license`_ for details.

.. _requirements: https://github.com/bartongroup/ProteoFAV/blob/master/requirements.txt
.. _license: https://github.com/bartongroup/ProteoFAV/blob/master/LICENSE.md
.. _issue tracker: https://github.com/bartongroup/ProteoFAV/issues
.. _docs: https://github.com/bartongroup/ProteoFAV/blob/master/docs/index.rst
.. _Pandas: http://pandas.pydata.org/
.. _Contributing Guide: https://github.com/bartongroup/ProteoFAV/wiki/Contributing-Guide
.. _Changelog: https://github.com/bartongroup/ProteoFAV/blob/master/CHANGELOG.rst
.. _Credits: https://github.com/bartongroup/ProteoFAV/blob/master/AUTHORS.rst
