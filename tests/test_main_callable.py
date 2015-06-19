#!/local/bin/python
# -*- coding: utf-8 -*-

'''
-------------------
custom_callables.py
-------------------

A module containing tests for the standard argparser custom callables I'm
using in :ref:`custom_callables.py <custom_callables_autodocs>`. These tests use
the standard `unittest <http://docs.python.org/2/library/unittest.html>`_
package and extend the :py:class:`unittest.TestCase` class.

.. moduleauthor:: Nick Schurch <nschurch@dundee.ac.uk>

:module_version: 1.0
:created_on: 2013-04-08

'''

__version__ = "1.0"

import unittest, os, tempfile, shutil, argparse, warnings
from utils import main_callable as cc


class TestParsingRoutines(unittest.TestCase):
    ''' Test the function defined in :ref:`custom_callables.py
    <custom_callables_autodocs>`.

    Each function is tested with existing and non-existent instances of the
    file or path type that it is designed for. We also test that providing non-
    string basic python objects produces string representations and doesn't
    case a failure.

    '''

    def setUp(self):
        ''' Initialize the framework for testing.

        Define and create a new system temporary directory with the
        `tempfile <http://docs.python.org/2/library/tempfile.html>`_ package
        for use as an existing directory. A file called *thisfileexists.txt*
        is created in this directory as an existing file. Two new paths that
        do not exist, one a fake sub-directory in the existing directory and
        the other a fake file in this fake directory, are constructed.
        '''

        # setup a new temp directory with a file in it
        self.existing_path = tempfile.mkdtemp()
        self.existing_file = "%s/thisfileexists.txt" % self.existing_path
        temp_file = open(self.existing_file, "w")
        temp_file.write("this exists")
        temp_file.close()

        # make up a mythical path and file based on the existing dir
        self.nonexisting_path = self.existing_path + "/testfornonexistingdir"
        self.nonexisting_file = self.nonexisting_path + "/thisdoesntexists.txt"

    def test_input_path(self):
        ''' Test the *input_path* function.

        =================    ==========================================
        Test object          Expectation
        =================    ==========================================
        Existing path        Return full path
        Non-existent path    Raise argparse.ArgumentTypeError exception
        Non-string arg       Raise TypeError exception
        =================    ==========================================
        '''

        result = cc.input_path(self.existing_path)
        msg = "Not finding existing path"
        self.assertEqual(result, self.existing_path, msg)

        # make sure this raises an argparse error
        self.assertRaises(argparse.ArgumentTypeError,
                          cc.input_path,
                          self.nonexisting_path)

        # make sure this raises an Type error if not given a string
        self.assertRaises(TypeError, cc.input_path, 5)
        self.assertRaises(TypeError, cc.input_path, [5])
        self.assertRaises(TypeError, cc.input_path, {"5": 5})

    def test_input_file(self):
        ''' Test the *input_file* function.

        =================    ==========================================
        Test object          Expectation
        =================    ==========================================
        Existing file        Return full file path
        Non-existent file    Raise argparse.ArgumentTypeError exception
        Non-string arg       Raise TypeError exception
        =================    ==========================================
        '''

        result = cc.input_file(self.existing_file)
        msg = "Not finding existing file"
        self.assertEqual(result, self.existing_file, msg)

        # make sure this raises an argparse error
        self.assertRaises(argparse.ArgumentTypeError,
                          cc.input_file,
                          self.nonexisting_file)

        # make sure this raises an Type error if not given a string
        self.assertRaises(TypeError, cc.input_file, 5)
        self.assertRaises(TypeError, cc.input_file, [5])
        self.assertRaises(TypeError, cc.input_file, {"5": 5})

    def test_output_path(self):
        ''' Test the *output_path* function.

        =================    ==========================================
        Test object          Expectation
        =================    ==========================================
        Existing path        Return full path
        Non-existent path    Create path & return full path
        Non-string arg       Raise TypeError exception
        =================    ==========================================
        '''

        result = cc.output_path(self.existing_path)
        msg = "Not finding existing path"
        self.assertEqual(result, self.existing_path, msg)

        # ignore expected warning messgaes
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = cc.output_path(self.nonexisting_path)

        msg = "Not creating the new path"
        self.assertEqual(result, self.nonexisting_path, msg)

        msg = "Not returning the correct path string"
        self.assertTrue(os.path.exists(self.nonexisting_path), msg)

        # make sure this raises an Type error if not given a string
        self.assertRaises(TypeError, cc.output_path, 5)
        self.assertRaises(TypeError, cc.output_path, [5])
        self.assertRaises(TypeError, cc.output_path, {"5": 5})

    def test_output_file(self):
        ''' Test the *output_file* function.

        =================    =============================================
        Test object          Expectation
        =================    =============================================
        Existing file        Return full file path
        Non-existent file    Create path (but not file) & return full path
        Non-string arg       Raise TypeError exception
        =================    =============================================
        '''

        result = cc.output_file(self.existing_file)
        msg = "Not finding existing file"
        self.assertEqual(result, self.existing_file, msg)

        # ignore expected warning messgaes
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = cc.output_path(self.nonexisting_file)

        msg = "Not returning the correct file string"
        self.assertEqual(result, self.nonexisting_file, msg)

        msg = "Not creating the new path"
        self.assertTrue(os.path.exists(self.nonexisting_path), msg)

        # make sure this raises an Type error if not given a string
        self.assertRaises(TypeError, cc.output_path, 5)
        self.assertRaises(TypeError, cc.output_path, [5])
        self.assertRaises(TypeError, cc.output_path, {"5": 5})

    def tearDown(self):
        ''' Remove testing framework.

        Deletes the created temporary directory and *thisfileexists.txt* file.
        '''

        # remove the tmpdir tree
        shutil.rmtree(self.existing_path)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestParsingRoutines)
    unittest.TextTestRunner(verbosity=2).run(suite)