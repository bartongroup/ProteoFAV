ProteoFAV
=========

**Protein Feature Aggregation and Variants**

|Pypi| |Build Status| |Documentation| |Python: versions| |License|

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
.. |License| image:: http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat
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

Getting Started
---------------

Dependencies
~~~~~~~~~~~~

ProteoFAV was developed to support Python 3.5+ and Pandas 0.20+.

Check `requirements`_ for specific requirements.

.. _requirements: https://github.com/bartongroup/ProteoFAV/blob/master/requirements.txt


Installation
~~~~~~~~~~~~

To install the stable release, run this command in your terminal:

.. code-block:: console

    $ pip install proteofav

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


Installing from source in a virtual environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Getting ProteoFAV:

.. code-block:: bash

    $ wget https://github.com/bartongroup/ProteoFAV/archive/master.zip -O ProteoFAV.zip
    $ unzip ProteoFAV.zip

    # alternatively, cloning the git repository
    $ git clone https://github.com/bartongroup/ProteoFAV.git


Installing With Conda:

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
    $ python path/to/ProteoFAV/setup.py install


Testing the installation
~~~~~~~~~~~~~~~~~~~~~~~~

Test dependencies should be resolved with:

.. code-block:: bash

    $ python path/to/ProteoFAV/setup.py develop --user


Run the Tests with:

.. code-block:: bash

    $ python path/to/ProteoFAV/setup.py test
    # or
    $ cd path/to/ProteoFAV/tests
    $ python -m unittest discover


ProteoFAV Configuration
~~~~~~~~~~~~~~~~~~~~~~~

ProteoFAV uses a configuration file `config.ini` where the user can specify the directory paths, as well as urls for commonly used data sources.

After installing run:

.. code-block:: bash

    $ proteofav-setup


Example Usage
-------------

Example usage is currently provided as a `Jupyter Notebook`, which can be viewed with the `GitHub's`_ file viewer or with the Jupyter `nbviewer`_.

You can download the Jupyter notebook from `GitHub`_ and test it with your ProteoFAV's installation.

.. _GitHub's: https://github.com/bartongroup/ProteoFAV/blob/master/Examples.ipynb
.. _nbviewer: https://nbviewer.jupyter.org/github/bartongroup/ProteoFAV/blob/master/Examples.ipynb
.. _GitHub: https://github.com/bartongroup/ProteoFAV


Contributing and Bug tracking
-----------------------------

Feel free to fork, clone, share and distribute. If you find any bugs or
issues please log them in the `issue tracker`_.

Before you submit your *Pull-requests* read the `Contributing Guide`_.

Credits
-------

See the `Credits`_


Changelog
---------

See the `Changelog`_


Licensing
---------

The MIT License (MIT). See `license`_ for details.

.. _requirements: https://github.com/bartongroup/ProteoFAV/blob/master/requirements.txt
.. _license: https://github.com/bartongroup/ProteoFAV/blob/master/LICENSE.md
.. _issue tracker: https://github.com/bartongroup/ProteoFAV/issues
.. _docs: https://github.com/bartongroup/ProteoFAV/blob/master/docs/index.rst
.. _Pandas: http://pandas.pydata.org/
.. _Contributing Guide: https://github.com/bartongroup/ProteoFAV/wiki/Contributing-Guide
.. _Changelog: https://github.com/bartongroup/ProteoFAV/blob/master/CHANGELOG.rst
.. _Credits: https://github.com/bartongroup/ProteoFAV/blob/master/AUTHORS.rst
