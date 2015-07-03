#!/local/bin/python
# -*- coding: utf-8 -*-

'''
--------------------
standard_parsers.py
--------------------

This module tests the :py:class:`argparse.ArgumentParser` instance returned by
:ref:`standard_parsers.py <standard_parsers_autodocs>` with a variety of
options and choices. These tests use the standard `unittest <http://docs.python
.org/2/library/unittest.html>`_ package and extend the
:py:class:`unittest.TestCase` class.

.. moduleauthor:: Nick Schurch <nschurch@dundee.ac.uk>

:module_version: 1.0
:created_on: 2013-04-08

'''

__version__ = "1.0"

import unittest, os, tempfile, shutil, warnings, argparse
import main_parser as sp


class TestStandardArgparse(unittest.TestCase):
    ''' Test the function defined in :ref:`standard_parsers.py
    <standard_parsers_autodocs>`.

    A standard :py:class:`argparse.ArgumentParser` instance is created from a
    call to :ref:`standard_parser <standard_parsers_autodocs>` and then is
    tested with different sets of arguments.

    '''

    def setUp(self):

        ''' Initialize the framework for testing.

        Creates a new parser and creates a new system temporary directory with
        the `tempfile <http://docs.python.org/2/library/tempfile.html>`_ package
        for use as an existing directory. A file called *thisfileexists.txt*
        is created in this directory as an existing file. Two new paths that
        do not exist, one a fake sub-directory in the existing directory and
        the other a fake file in this fake directory, are constructed.
        '''

        # set up an existing path and file
        self.existing_path = tempfile.mkdtemp()
        self.existing_file = "%s/thisfileexists.txt" % self.existing_path
        temp_file = open(self.existing_file, "w")
        temp_file.write("this exists")
        temp_file.close()

        # make up a mythical path and file based on the existing dir
        self.nonexisting_path = self.existing_path + "/testfornonexistingdir"
        self.nonexisting_file = self.nonexisting_path + "/thisdoesntexists.txt"

    def tearDown(self):

        ''' Remove testing framework.

        Deletes the created temporary directory and *thisfileexists.txt* file.
        '''

        shutil.rmtree(self.existing_path)

    def test_standard_parser_returns(self):

        ''' Test that the standard_parser function returns the correct objects

        =================    ================================================
        Test object          Expectation
        =================    ================================================
        valid arguments      Returns a :py:class:`argparse.ArgumentParser`
                             instance & 2 lists; 1) list of positional
                             arguments & defaults (should be None), 2)list of
                             keyword arguments & defaults.
        =================    ================================================
        '''

        msg = "standard parser not returning the correct number of things"
        parser_tup = sp.standard_parser(__version__)
        self.assertTrue(len(parser_tup) == 3, msg)

        opts, pos_args, kw_args = parser_tup

        msg = "standard parser not returning the correct optional argument list"
        # test kw_args
        self.assertTrue(len(kw_args) == 3, msg)

        foundVerbose = False
        foundLog = False
        foundTmpdir = False

        for arg in kw_args:
            if arg[0] == "verbose":
                foundVerbose = True
                self.assertTrue(arg == ("verbose", False), msg)
            elif arg[0] == "log":
                foundLog = True
                self.assertTrue(arg == ("log", None), msg)
            elif arg[0] == "tmpdir":
                foundTmpdir = True

        self.assertTrue(foundVerbose, msg)
        self.assertTrue(foundLog, msg)
        self.assertTrue(foundTmpdir, msg)

        # test pos_args
        msg = "standard parser not returning the correct positional " \
              "argument list"
        self.assertTrue(len(pos_args) == 2, msg)
        self.assertTrue(("infile", None) in pos_args, msg)
        self.assertTrue(("outfile", None) in pos_args, msg)

        # test opts
        msg = "standard parser not an argparse.ArgumentParser instance"
        self.assertTrue(isinstance(opts, argparse.ArgumentParser), msg)

    def test_standard_parser_infile(self):

        ''' Test the parser with existing and non-existent infile arguments

        =================    ==========================
        Test object          Expectation
        =================    ==========================
        Existing file        Stores the filename
        Non-existent file    Raise SystemExit exception
        =================    ==========================
        '''

        parser, pos_args, kw_args = sp.standard_parser(__version__)

        # test infile with existing file
        fakeargs = [self.existing_file,
                    self.existing_file,
                    "-l", self.existing_file
                    ]

        msg = "existing infile not being handled correctly"
        opts = parser.parse_args(fakeargs)
        self.assertEqual(opts.infile, self.existing_file, msg)

        # test infile with non-existent file
        fakeargs = [self.nonexisting_file,
                    self.existing_file,
                    "-l", self.nonexisting_file
                    ]

        self.assertRaises(SystemExit,
                          parser.parse_args,
                          fakeargs)

    def test_standard_parser_outfile(self):

        ''' Test the parser with existing and non-existent outfile arguments

        =================    ========================================
        Test object          Expectation
        =================    ========================================
        Existing file        Stores the filename
        Non-existent file    Creates the path and stores the filename
        =================    ========================================
        '''

        # test outfile with existing file
        parser, pos_args, kw_args = sp.standard_parser(__version__)
        fakeargs = [self.existing_file,
                    self.existing_file,
                    "-l", self.existing_file
                    ]

        msg = "existing outfile not being handled correctly"
        opts = parser.parse_args(fakeargs)
        self.assertEqual(opts.outfile, self.existing_file, msg)

        # test outfile with non-existent file
        fakeargs = [self.existing_file,
                    self.nonexisting_file,
                    "-l", self.existing_file
                    ]

        msg = "non-existent outfile and dir not being handled correctly"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            opts = parser.parse_args(fakeargs)
            self.assertTrue(os.path.exists(self.nonexisting_path), msg)
            self.assertEqual(opts.outfile, self.nonexisting_file, msg)

    def test_standard_parser_log(self):

        ''' Test the parser with existing and non-existent -l|--log arguments

        =================    ========================================
        Test object          Expectation
        =================    ========================================
        Existing file        Stores the filename
        Non-existent file    Creates the path and stores the filename
        =================    ========================================
        '''

        # test logfile with existing file
        parser, pos_args, kw_args = sp.standard_parser(__version__)
        fakeargs = [self.existing_file,
                    self.existing_file,
                    "-l", self.existing_file
                    ]

        msg = "existing logfile not being handled correctly"
        opts = parser.parse_args(fakeargs)
        self.assertEqual(opts.log, self.existing_file, msg)

        # test logfile with non-existant dir
        fakeargs = [self.existing_file,
                    self.existing_file,
                    "-l", self.nonexisting_file
                    ]

        msg = "non-existant log file and dir not being handled correctly"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            opts = parser.parse_args(fakeargs)
            self.assertTrue(os.path.exists(self.nonexisting_path), msg)
            self.assertEqual(opts.log, self.nonexisting_file, msg)

    def test_standard_parser_verbosity(self):

        ''' Test the parser with a -v|--verbose argument

        =================    ========================================
        Test object          Expectation
        =================    ========================================
        no -v                Stores False
        -v                   Stores True
        =================    ========================================
        '''

        # test verbosity
        parser, pos_args, kw_args = sp.standard_parser(__version__)
        fakeargs = [self.existing_file,
                    self.existing_file,
                    "-l", self.existing_file
                    ]

        msg = "verbosity not being parsed correctly"
        opts = parser.parse_args(fakeargs)
        self.assertFalse(opts.verbose, msg)

        fakeargs = [self.existing_file,
                    self.existing_file,
                    "-l", self.existing_file,
                    "-v"
                    ]

        opts = parser.parse_args(fakeargs)
        self.assertTrue(opts.verbose, msg)

    def test_standard_parser_tmpdir(self):

        ''' Test the parser with existing and non-existent --tmpdir arguments

        ==============================  ========================================
        Test object                     Expectation
        ==============================  ========================================
        no --tmpdir                     Create & store a system defined temp dir
        --tmpdir                        Create & store a system defined temp dir
        --tmpdir *<existing path>*      Store *<path>* as the temp dir
        --tmpdir *<non-existent path>*  Create & store *<path>* as the temp dir
        ==============================  ========================================
        '''

        parser, pos_args, kw_args = sp.standard_parser(__version__)
        fakeargs = [self.existing_file,
                    self.existing_file,
                    "-l", self.existing_file
                    ]

        msg = "omitted tmpdir arg not being handled correctly"
        opts = parser.parse_args(fakeargs)
        self.assertTrue(os.path.exists(opts.tmpdir), msg)

        fakeargs = [self.existing_file,
                    self.existing_file,
                    "-l", self.existing_file,
                    "--tmpdir"
                    ]

        msg = "unspecified tmpdir arg not being handled correctly"
        opts = parser.parse_args(fakeargs)
        self.assertTrue(os.path.exists(opts.tmpdir), msg)

        # test tmpdir with existing but non-random dir
        fakeargs = [self.existing_file,
                    self.existing_file,
                    "-l", self.existing_file,
                    "--tmpdir", self.existing_path
                    ]

        msg = "existing tmpdir not being handled correctly"
        opts = parser.parse_args(fakeargs)
        self.assertEqual(opts.tmpdir, self.existing_path, msg)

        # test tmpdir with non-existant dir
        fakeargs = [self.existing_file,
                    self.existing_file,
                    "-l", self.existing_file,
                    "--tmpdir", self.nonexisting_path
                    ]

        msg = "non-existent tmpdir not being handled correctly"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            opts = parser.parse_args(fakeargs)
            self.assertEqual(opts.tmpdir, self.nonexisting_path, msg)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestStandardArgparse)
    unittest.TextTestRunner(verbosity=2).run(suite)
