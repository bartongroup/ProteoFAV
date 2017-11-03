.. highlight:: shell

=============================
Contributing and Bug tracking
=============================

Feel free to fork, clone, share and distribute. If you find any bugs or issues please log them in the `issue tracker`_.

Before you submit your Pull-requests read the `Contributing Guide`_ below.

.. _issue tracker: https://github.com/bartongroup/ProteoFAV/issues
.. _Contributing Guide: contributing.html#contributing-guide


Contributing Guide
------------------

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given. You can
contribute in many ways.

Types of Contributions
----------------------

Please check GitHub issue tracker for:

- bugs
- features
- enhancements

Style
-----

Wrap line over 99 (79 if possible). Lines with more than 79 characters
are harder to read in the command line. Else, we try to keep the
pep8 standards. Private functions start with an underline.

Naming Conventions
------------------

Functions that: get data from a resources should start with **fetch/\_fetch**.
To load and process python objects **select/\_select**.

Python code is also documentation, hence explicit variable names are
recommended: **ensembl\_ptn\_id** instead **identifier**. Generally,
for a function input parameter **identifier** can be used if the docstring
clearly defines which type of Accession Identifier should be used.

If the resource has its own name for a field or value, try to keep, for
consistency.

Docstrings
----------

Use the following template:

::

    def df_encoder(data, descriptor=None, col_names=None):
        """Encode a pandas DataFrame with a descriptor. Similar to one hot
            encoding, however it preserves encoding between columns.

        :param descriptor: dict like descriptor to be applied to the columns
        :type descriptor:
        :param data: pandas DataFrame with categorical data
        :type data: pandas.DataFrame
        :param col_names: names for the new columns
        :type col_names: list of [str,]
        :return: table with encoded data
        :rtype : pd.DataFrame

        :Example:

            >>> import pandas as pd
            >>> df = pd.DataFrame(map(list, ['ABC', 'AAA', 'ACD']))
            >>> print(df_encoder(df))
               0_A  1_A  2_A  0_B  1_B  2_B  0_C  1_C  2_C  0_D  1_D  2_D
            0    1    0    0    0    0    1    0    0    0    0    1    0
            1    1    0    0    0    1    0    0    0    1    0    0    0
            2    1    0    0    0    0    0    1    0    0    0    0    1

          .. note:: Use an external descriptor to make sure you descriptor is
          replicated
         """

If the function returns a pandas.DataFrame, is good practice to add
    which columns: dtype you expect, so we can keep track of it. Since
    there is no common standard for Pandas column annotation, we are
    still deciding the best approach:

::

    def fetch_ensembl_variants(ensembl_ptn_id, feature=None):
    """Queries the Ensembl API for germline variants (mostly dbSNP) and somatic
    (mostly COSMIC) based on Ensembl Protein identifiers (e.g. ENSP00000326864).

    :param ensembl_ptn_id: Ensembl acession to a protein: ENSP00000XXXXXX
    :return: table[Parent: str,
                   allele: str,
                   clinical_significance: list,
                   codons: str,
                   end: int,
                   feature_type: str,
                   id: str,
                   minor_allele_frequency: float ,
                   polyphen: float,
                   residues: str,
                   seq_region_name: str,
                   sift: float,
                   start: int,
                   translation: str,
                   type: str ]
    :rtype: pandas.DataFrame
    """

Columns type pragma
-------------------

Culumn type normalisation is a central issue in ProteoFAV. There is no simple
way to make column type consintent across all data files. Some pragmatic
rules to deal with NANs (Not an number) in non float columns are defined
here, but open to change. NAN in Python are always floats. If one has to
operate with integers or string, it must eliminate the NAN’s, and in
ProteoFAV we use the following rules:

* If is a sequence index: -9999
* If is a sequence column NAN’s: ‘X’
* If is another string column: '' (empty string)

Testing
-------

Doctests are not mandatory, but tests are. Tests are located in `/tests`
and we use standard Unittest setup.
