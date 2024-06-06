# -*- coding: utf-8 -*-
#
#
# TheVirtualBrain-Framework Package. This package holds all Data Management, and 
# Web-UI helpful to run brain-simulations. To use it, you also need to download
# TheVirtualBrain-Scientific Package (for simulators). See content of the
# documentation-folder for more details. See also http://www.thevirtualbrain.org
#
# (c) 2012-2024, Baycrest Centre for Geriatric Care ("Baycrest") and others
#
# This program is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this
# program.  If not, see <http://www.gnu.org/licenses/>.
#
#
#   CITATION:
# When using The Virtual Brain for scientific publications, please cite it as explained here:
# https://www.thevirtualbrain.org/tvb/zwei/neuroscience-publications
#
#
"""
.. moduleauthor:: Bogdan Neacsa <bogdan.neacsa@codemart.ro>
"""

import os
import tvb_data

from tvb.adapters.datatypes.db.connectivity import ConnectivityIndex
from tvb.tests.framework.core.base_testcase import TransactionalTestCase
from tvb.adapters.visualizers.connectivity import ConnectivityViewer
from tvb.tests.framework.core.factory import TestFactory


class TestConnectivityViewer(TransactionalTestCase):
    """
    Unit-tests for Connectivity Viewer.
    """

    def transactional_setup_method(self):
        """
        Sets up the environment for running the tests;
        creates a test user, a test project, a connectivity and a surface;
        imports a CFF data-set
        """

        self.test_user = TestFactory.create_user("UserCVV")
        self.test_project = TestFactory.create_project(self.test_user)

        zip_path = os.path.join(os.path.dirname(tvb_data.__file__), 'connectivity', 'connectivity_66.zip')
        TestFactory.import_zip_connectivity(self.test_user, self.test_project, zip_path)
        self.connectivity_index = TestFactory.get_entity(self.test_project, ConnectivityIndex)
        assert self.connectivity_index is not None

    def test_launch(self):
        """
        Check that all required keys are present in output from BrainViewer launch.
        """
        viewer = ConnectivityViewer()
        view_model = viewer.get_view_model_class()()
        view_model.connectivity = self.connectivity_index.gid
        result = viewer.launch(view_model)
        expected_keys = ['weightsMin', 'weightsMax', 'urlWeights', 'urlVertices',
                         'urlTriangles', 'urlTracts', 'urlPositions', 'urlNormals',
                         'rightHemisphereJson', 'raysArray', 'rayMin', 'rayMax', 'positions',
                         'leftHemisphereJson', 'connectivity_entity', 'bothHemisphereJson']
        for key in expected_keys:
            assert key in result
