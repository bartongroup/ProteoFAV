ProteoFAV's Documentation
=========================

ProteoFAV (Protein Feature Aggregation and Variants) is an open-source framework for simple and fast integration of protein structure data with sequence annotations and genetic variation.

ProteoFAV is a Python module that address the challenge of cross-mapping protein structures and protein sequences,
allowing for protein structures to be annotated with sequence features. It implements methods for working with
protein structures (via mmCIF, PDB, PDB Validation, DSSP and SIFTS files), sequence Features (via UniProt GFF annotations) and
genetic variants (via UniProt/EBI Proteins API and Ensembl REST API). Cross-mapping of structure and sequence is
performed with the aid of SIFTS.

ProteFAV relies heavily in the `Pandas`_ library to quickly load data into DataFrames for fast
data exploration and analysis. Structure and sequence
data are parsed/fetched onto Pandas DataFrames that are then merged-together (collapsed) onto a
single DataFrame.

.. _Pandas: http://pandas.pydata.org/


Table of Contents
=================

.. toctree::
   :maxdepth: 5

   getting_started
   example_usage
   contributing
   authors
   changelog
   proteofav_docs


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
