# coding=utf-8
"""DockWidget test.

.. note:: This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 2 of the License, or
     (at your option) any later version.

"""

__author__ = 'miguel.canon@tnc.org'
__date__ = '2024-09-23'
__copyright__ = 'Copyright 2024, The Nature Conservancy'

import unittest

from qgis.PyQt.QtGui import QDockWidget

from Siga_Step_01_Config_dockwidget import Siga_Step_01_ConfigDockWidget

from utilities import get_qgis_app

QGIS_APP = get_qgis_app()


class Siga_Step_01_ConfigDockWidgetTest(unittest.TestCase):
    """Test dockwidget works."""

    def setUp(self):
        """Runs before each test."""
        self.dockwidget = Siga_Step_01_ConfigDockWidget(None)

    def tearDown(self):
        """Runs after each test."""
        self.dockwidget = None

    def test_dockwidget_ok(self):
        """Test we can click OK."""
        pass

if __name__ == "__main__":
    suite = unittest.makeSuite(Siga_Step_01_ConfigDialogTest)
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)

