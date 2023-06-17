
# -*- coding: utf-8 -*-
#
#
# TheVirtualBrain-Scientific Package. This package holds all simulators, and
# analysers necessary to run brain-simulations. You can use it stand alone or
# in conjunction with TheVirtualBrain-Framework Package. See content of the
# documentation-folder for more details. See also http://www.thevirtualbrain.org
#
# (c) 2012-2023, Baycrest Centre for Geriatric Care ("Baycrest") and others
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
.. moduleauthor:: Abhijit Deo <f20190041@goa.bits-pilani.ac.in>
"""


#TODO : change the record id to the latest when done testing on local machine. :) :wq

from tvb.datasets import Zenodo, Record
from pathlib import Path
from tvb.tests.library.base_testcase import BaseTestCase


class TestZenodo(BaseTestCase):

    def test_get_record(self):

        zenodo = Zenodo()
        rec = zenodo.get_record("4263723")

        assert type(rec) == Record 
        assert rec.data["doi"] == "10.5281/zenodo.4263723"

        del rec 
        del zenodo


    def test_get_versions(self):

        zenodo = Zenodo()
        versions = zenodo.get_versions_info("3491055")

        assert type(versions) == dict
        assert versions == {'2.0.1': '3497545', '1.5.9.b': '3474071', '2.0.0': '3491055', '2.0.3': '4263723', '2.0.2': '3688773', '1.5.9': '3417207', '2.7': '7574266'}

        del zenodo
        del versions

class TestRecord(BaseTestCase):


    def test_download(self):

        zen = Zenodo()

        rec = zen.get_record("4263723")


        rec.download()
        print(rec.file_loc)
        for file_name, file_path in rec.file_loc.items():
            assert Path(file_path).is_file()


