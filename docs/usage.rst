=====
Usage
=====

ProteoFAV works both as a module and a command line tool.

From the command line use::

    python PycharmProjects/ProteoFAV/ --pdb=2w4o --chain=A test.csv


To use ProteoFAV in a project::

    import proteofav

To import a protected method - do it at your own risk - use::

    from proteofav.uniprot import _uniprot_info
