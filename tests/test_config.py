# -*- coding: utf-8 -*-

import os
import sys
import logging
import unittest

try:
    from mock import patch
except ImportError:
    from unittest.mock import patch

from proteofav.config import Defaults
from proteofav.config import defaults as config

mock_config = """\
[Global]
db_pdb = /pdb_dir/

[Other]
test = /test/value
"""


@patch("proteofav.config.defaults.db_pdb", "/pdb/")
@patch("proteofav.config.defaults.db_mmcif", "/mmcif/")
class TestConfig(unittest.TestCase):
    """Test the Config methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.config = config
        self.defaults = Defaults

    def tearDown(self):
        """Remove testing framework."""

        self.config = None
        self.defaults = None

    def test_loading_config_defaults(self):
        config = self.config
        self.assertTrue(hasattr(config, 'db_pdb'))
        self.assertEqual(config.db_pdb, "/pdb/")
        self.assertTrue(hasattr(config, 'db_mmcif'))
        self.assertEqual(config.db_mmcif, "/mmcif/")
        self.assertFalse(hasattr(config, 'test'))

    def test_updating_config_defaults(self):
        config = self.config
        config.db_pdb = '/new_value/'
        self.assertNotEqual(config.db_pdb, "/pdb/")
        self.assertEqual(config.db_pdb, "/new_value/")

    def test_deleting_entry_config(self):
        config = self.config
        del config.db_pdb
        self.assertFalse(hasattr(config, 'db_pdb'))

    def test_loading_config_from_file(self):
        # using mocking config in this case
        self.config.db_pdb = '/tmp/'
        new_config_file = os.path.join(self.config.db_pdb,
                                       "mock_config.ini")
        with open(new_config_file, 'w') as out:
            out.write(mock_config)
        config = self.defaults(config_file=new_config_file)
        self.assertFalse(hasattr(config, 'some_name'))
        self.assertTrue(hasattr(config, 'test'))
        self.assertEqual(config.test, "/test/value")
        os.remove(new_config_file)

    def test_update_config_with_new_file(self):
        # load initial config
        config = self.defaults()
        self.assertFalse(hasattr(config, 'some_name'))
        self.assertFalse(hasattr(config, 'test'))

        # write mock config to file
        self.config.db_pdb = '/tmp/'
        new_config_file = os.path.join(self.config.db_pdb,
                                       "mock_config.ini")
        with open(new_config_file, 'w') as out:
            out.write(mock_config)

        # update config by loading the new config
        config.update(config_file=new_config_file)
        self.assertFalse(hasattr(config, 'some_name'))
        self.assertTrue(hasattr(config, 'test'))
        self.assertEqual(config.test, "/test/value")
        os.remove(new_config_file)

    def test_write_default_config(self):
        config = self.defaults()
        # write new config to file
        config.db_pdb = '/tmp/'
        new_config_file = os.path.join(config.db_pdb, 'tmp.config')
        config.write(new_config_file)

        # re-load default config based on the file
        config = self.defaults(config_file=new_config_file)
        self.assertTrue(hasattr(config, 'db_pdb'))
        self.assertTrue(hasattr(config, 'db_mmcif'))
        self.assertFalse(hasattr(config, 'test'))
        os.remove(new_config_file)


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("proteofav.config").setLevel(logging.CRITICAL)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestConfig)
    unittest.TextTestRunner(verbosity=2).run(suite)
