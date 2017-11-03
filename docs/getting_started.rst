.. highlight:: shell

===============
Getting Started
===============

Dependencies
------------

ProteoFAV was developed to support Python 3.5+ and Pandas 0.20+.

Check `requirements`_ for specific requirements.

.. _requirements: https://github.com/bartongroup/ProteoFAV/blob/master/requirements.txt


Installation
------------

To install the stable release, run this command in your terminal:

.. code-block:: console

    $ pip install proteofav

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


Installing from source in a virtual environment
-----------------------------------------------

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
------------------------

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
-----------------------

ProteoFAV uses a configuration file `config.ini` where the user can specify the directory paths, as well as urls for commonly used data sources.

After installing run:

.. code-block:: bash

    $ proteofav-setup
