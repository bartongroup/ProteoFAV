=========
Changelog
=========

0.2 beta (Nov-2017)
-------------------

* New interfaces to Data (via class methods 'read', 'write', 'download', 'select' and 'fetch')
* New modular approach for merging Pandas Tables (via methods 'merge' and 'generate')
* Parsing of PDB-formatted files and writing in both mmCIF and PDB formats
* Improved aggregation and flattening (normalisation) of genetic variants from UniProt Proteins API and Ensembl REST API
* Updated Documentation and CI


0.1 alpha (2015-2017)
---------------------

* Parsing methods for various file formats (mmCIF, DSSP, SIFTS, PDB Validation and GFF)
* Fetching methods accessing data from various APIs (PDBe API, Ensembl REST API, UniProt API, etc.)
* A main `merge_tables` function and CLI for generating a 'big' Pandas DataFrame with aggregated Structural Data, Sequence Features and Variants


0.0 (29-05-2015)
----------------

* ProteoFAV was created.
