# -*- coding: utf-8 -*-
#
#
# TheVirtualBrain-Framework Package. This package holds all Data Management, and
# Web-UI helpful to run brain-simulations. To use it, you also need do download
# TheVirtualBrain-Scientific Package (for simulators). See content of the
# documentation-folder for more details. See also http://www.thevirtualbrain.org
#
# (c) 2012-2022, Baycrest Centre for Geriatric Care ("Baycrest") and others
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
# When using The Virtual Brain for scientific publications, please cite it as follows:
#
#   Paula Sanz Leon, Stuart A. Knock, M. Marmaduke Woodman, Lia Domide,
#   Jochen Mersmann, Anthony R. McIntosh, Viktor Jirsa (2013)
#       The Virtual Brain: a simulator of primate brain network dynamics.
#   Frontiers in Neuroinformatics (7:10. doi: 10.3389/fninf.2013.00010)
#
#

"""
Classes and enums for the 3 analysis levels relevant for the image processing pipeline.

.. moduleauthor:: Robert Vincze <robert.vincze@codemart.ro>
"""

from tvb.basic.neotraits.api import HasTraits, Attr, NArray, Range, TVBEnum
from tvb.core.neotraits.traits_with_parameters import TraitsWithParameters


class ParcellationOptionsEnum(str, TVBEnum):
    AAL_PARC = "aal"
    AAL2_PARC = "aal2"
    BRAINNETOME_PARC = "brainnetome246fs"
    CRADDOCK200_PARC = "craddock200"
    CRADDOCK400_PARC = "craddock400"
    DESIKAN_PARC = "desikan"
    DESTRIEUX_PARC = "destrieux"
    HCPMMP1_PARC = "hcpmmp1"
    PERRY512_PARC = "perry512"
    YEO7fs_PARC = "yeo7fs"
    YEO7mni_PARC = "yeo7mni"
    YEO17fs_PARC = "yeo17fs"
    YEO17mni_PARC = "yeo17mni"


class TemplateRegOptionsEnum(str, TVBEnum):
    ANTS_TEMPLATE_REG = "ants"
    FSL_TEMPLATE_REG = "fsl"


class PipelineAnalysisLevel(TraitsWithParameters):
    pass


class PreprocAnalysisLevel(PipelineAnalysisLevel):
    parameters = Attr(
        field_type=dict,
        label="Linear Parameters",
        default=lambda: {"t1w_preproc_path": ""})

    def __str__(self):
        return "Preproc Analysis Level"


class ParticipantAnalysisLevel(PipelineAnalysisLevel):
    parameters = Attr(
        field_type=dict,
        label="Linear Parameters",
        default=lambda: {"t1w_preproc_path": "", "parcellation": ParcellationOptionsEnum.AAL_PARC, "stream_lines": 1,
                         "template_reg": TemplateRegOptionsEnum.ANTS_TEMPLATE_REG})

    def __str__(self):
        return "Participant Analysis Level"


class GroupAnalysisLevel(PipelineAnalysisLevel):
    parameters = Attr(
        field_type=dict,
        label="Linear Parameters",
        default=lambda: {"group_participant_label": "default", "session_label": "default"})

    def __str__(self):
        return "Group Analysis Level"
