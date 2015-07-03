#!/local/bin/python
# -*- coding: utf-8 -*-

'''
-------------------
standard_logging.py
-------------------

This module tests the :py:class:`logging.Logger` instance returned by
:ref:`standard_logging.py <standard_logging_autodocs>` with a variety of
options and choices. These tests use the standard `unittest <http://docs.python
.org/2/library/unittest.html>`_ package and extend the
:py:class:`unittest.TestCase` class.

.. moduleauthor:: Nick Schurch <nschurch@dundee.ac.uk>

:module_version: 1.0
:created_on: 2013-04-08

'''

__version__ = "1.0"

import unittest, os, tempfile, shutil, warnings, argparse, re, logging, sys
import main_logger as sl
import main_parser as sp


class TestStandardLogging(unittest.TestCase):
    '''Test the functions defined in :ref:`standard_logging.py
    <standard_logging_autodocs>`.

    A standard :py:class:`logging.Logger` instance is created from a
    call to :ref:`standard_logging <standard_logging_autodocs>` and then is
    tested with different sets of arguments.

    '''

    def setUp(self):
        '''Initialize the framework for testing.

        Creates a new system temporary directory with the `tempfile
        <http://docs.python.org/2/library/tempfile.html>`_ package
        as a place to generate our test log file, and generate a name for the
        test log file. Also generates :py:class:`argparse.ArgumentParser`
        instance from :ref:`standard_parsers.py <standard_parsers_autodocs>`,
        a fake command line and arguments.

        '''

        self.existing_path = tempfile.mkdtemp()
        self.fake_infile = "%s/thisfileexists.txt" % self.existing_path
        temp_file = open(self.fake_infile, "w")
        temp_file.write("this exists")
        temp_file.close()
        self.fake_outfile = "%s/thisfileexists.txt" % self.existing_path
        self.log_file = "%s/test.log" % self.existing_path

        parser_tup = sp.standard_parser(__version__)
        self.parser, self.pos_args, self.kw_args = parser_tup

        self.parser_args = [self.fake_infile,
                            self.fake_outfile,
                            "-l", self.log_file,
                            "-v"
                            ]

        self.args = self.parser.parse_args(self.parser_args)
        self.fake_cmd_line = sys.argv + self.parser_args

    def tearDown(self):
        '''Remove testing framework.

        Deletes the created temporary directory and *thisfileexists.txt* file.
        '''
        logging.shutdown()
        # shutil.rmtree(self.existing_path)

    def test_standard_logger_attributes(self):
        '''Test for graceful fail when not given required arguments

        =====================    ==========================
        Test object              Expectation
        =====================    ==========================
        no *verbose* argument    Raise AttributeError
        no *log* argument        Raise AttributeError
        =====================    ==========================
        '''

        parser = argparse.ArgumentParser()
        parser.add_argument('--log')
        args = parser.parse_args(["--log", self.log_file])

        self.assertRaises(AttributeError,
                          sl.standard_logger,
                          __version__, sys.argv, args, [], [("log", None)]
                          )

        parser = argparse.ArgumentParser()
        parser.add_argument('-verbose', action="store_true")
        args = parser.parse_args(["-verbose"])

        self.assertRaises(AttributeError,
                          sl.standard_logger,
                          __version__, sys.argv, args, [], [("log", None)]
                          )

    def test_standard_logger_basic_logging(self):
        '''Test for open logfile & that script command-line is correctly written

        ============================  =========================================
        Test object                   Expectation
        ============================  =========================================
        valid command-line arguments  Create file and write welcome message,
                                      time date and calling command line to it.
        ============================  =========================================
        '''

        self.logger = sl.standard_logger(__version__,
                                         self.fake_cmd_line,
                                         self.args,
                                         self.pos_args,
                                         self.kw_args)

        msg = "log file is not being created correctly"
        self.assertTrue(os.path.isfile(self.log_file), msg)

        logfile = open(self.log_file)
        log_data = logfile.readlines()
        logfile.close()

        msg = "command-line call not being registered to the log file"
        match_critereon = "^.+%s$" % " ".join(self.fake_cmd_line)
        self.assertTrue(any(re.match(match_critereon, line)
                            for line in log_data
                            ),
                        msg)

    def test_standard_logger_parsed_arguments(self):
        '''Test for logging of parsed arguments

        ============================  =========================================
        Test object                   Expectation
        ============================  =========================================
        valid command-line arguments  Log the full set of script arguments
                                      parsed byt he parser. Include defaults
                                      and make sure they are marked.
        ============================  =========================================
        '''

        self.logger = sl.standard_logger(__version__,
                                         self.fake_cmd_line,
                                         self.args,
                                         self.pos_args,
                                         self.kw_args)

        logfile = open(self.log_file)
        log_data = logfile.readlines()
        logfile.close()

        msg = "parsed infile argument is not being registered to the log " \
              "file correctly"
        infile_match = "^infile.+:  %s.*$" % self.fake_infile
        self.assertTrue(any(re.match(infile_match, line.strip())
                            for line in log_data
                            ),
                        msg)

        msg = "parsed outfile argument is not being registered to the log " \
              "file correctly"
        outfile_match = "^outfile.+:  %s$" % self.fake_outfile
        self.assertTrue(any(re.match(outfile_match, line.strip())
                            for line in log_data
                            ),
                        msg)

        msg = "parsed log argument is not being registered to the log " \
              "file correctly"
        logfile_match = "^--log.+:  %s$" % self.log_file
        self.assertTrue(any(re.match(logfile_match, line.strip())
                            for line in log_data
                            ),
                        msg)

        msg = "parsed verbose argument is not being registered to the log " \
              "file correctly"
        verbose_match = "^--verbose.+:  .+$"
        self.assertTrue(any(re.match(verbose_match, line.strip())
                            for line in log_data
                            ),
                        msg)

        print self.log_file
        msg = "parsed default tmpdir argument is not being registered to the " \
              "log file correctly"
        tmpdir_match = "^--tmpdir.+:.+\(default\)$"
        self.assertTrue(any(re.match(tmpdir_match, line.strip())
                            for line in log_data
                            ),
                        msg)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestStandardLogging)
    unittest.TextTestRunner(verbosity=2).run(suite)
