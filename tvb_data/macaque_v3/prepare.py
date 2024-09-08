# -*- coding: utf-8 -*-
#
#
# TheVirtualBrain-Framework Package. This package holds all Data Management, and 
# Web-UI helpful to run brain-simulations. To use it, you also need do download
# TheVirtualBrain-Scientific Package (for simulators). See content of the
# documentation-folder for more details. See also http://www.thevirtualbrain.org
#
# (c) 2012-2017, Baycrest Centre for Geriatric Care ("Baycrest") and others
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

import os
import numpy
import nibabel as nib
from zipfile import ZipFile


def alter_tracts_and_weights():
    tract_lengths = numpy.loadtxt("Connectivity/tract_lengths_orig.txt")
    weights = numpy.loadtxt("Connectivity/weights_orig.txt")

    # add row of 0's at index 83
    tract_lengths = numpy.vstack((tract_lengths, numpy.zeros((1, 82))))
    weights = numpy.vstack((weights, numpy.zeros((1, 82))))

    # add row of 0's at index 41
    tract_lengths = numpy.vstack((tract_lengths[:41, ], numpy.zeros((1, 82)), tract_lengths[41:, ]))
    weights = numpy.vstack((weights[:41, ], numpy.zeros((1, 82)), weights[41:, ]))

    # add column of 0's at index 83
    tract_lengths = numpy.hstack((tract_lengths, numpy.zeros((84, 1))))
    weights = numpy.hstack((weights, numpy.zeros((84, 1))))

    # add column of 0's at index 41
    tract_lengths = numpy.hstack((tract_lengths[:, :41], numpy.zeros((84, 1)), tract_lengths[:, 41:]))
    weights = numpy.hstack((weights[:, :41], numpy.zeros((84, 1)), weights[:, 41:]))

    return tract_lengths, weights


def prepare_rm():
    original_rm = numpy.loadtxt("RM/surf_labels.txt")
    # We Manually unk_l and unk_R into positions 41 and 83 in the TVBmacaque_RM_LUT file received from Kelly
    labels_map = numpy.loadtxt("RM/TVBmacaque_RM_LUT.txt", skiprows=1, usecols=[0], dtype=numpy.int32)
    labels_text = numpy.loadtxt("RM/TVBmacaque_RM_LUT.txt", skiprows=1, usecols=[1], dtype=numpy.str_)
    vertices_x = numpy.loadtxt("Surface/vertices.txt", usecols=[0], dtype=numpy.float64)

    reverted_map = dict()
    for i in range(labels_map.size):
        reverted_map[labels_map[i]] = i
    print(reverted_map)

    n_surface = original_rm.size
    print(n_surface, original_rm.shape)
    print(numpy.count_nonzero(original_rm < 2))

    # values that don't belong to any region will be assigned to the unknown regions
    final_rm = []
    for i in range(n_surface):
        int(original_rm[i])
        if int(original_rm[i]) in reverted_map:
            final_rm.append(reverted_map[int(original_rm[i])])
        else:
            # if the Y coordinate of the vertex is positive, we will assign it to the unknown right region,
            # otherwise to the left one
            # print(i, vertices_x[i], original_rm[i])
            if (vertices_x[i] > 0):
                final_rm.append(41)
            else:
                final_rm.append(83)

    print(numpy.min(final_rm), numpy.max(final_rm))
    final_rm = numpy.array(final_rm, dtype=numpy.int32)
    numpy.savetxt("regionMapping_147k_84.txt", final_rm)

    maps = []
    for i in range(labels_map.size):
        maps.append((labels_map[i], labels_text[i], i))
    maps = numpy.array(maps, dtype=numpy.str_)
    numpy.savetxt("mapping_FS_84.txt", maps, delimiter=" ", fmt="%s")

    # alter and save weights and tract_lengths adjacency matrices
    tract_lengths, weights = alter_tracts_and_weights()
    numpy.savetxt("Connectivity/tract_lengths.txt", tract_lengths, delimiter=" ")
    numpy.savetxt("Connectivity/weights.txt", weights, delimiter=" ")


def prepare_centres():
    with open("RM/TVBmacaque_RM_LUT.txt", "r") as f:
        labels_text = [x.split() for x in f.readlines()]
    centers_ordered = [row[-1] for row in labels_text[1:]]

    centre_keys = numpy.loadtxt("Connectivity/centres_orig.txt", usecols=[0], dtype=numpy.str_)
    centre_positions = numpy.loadtxt("Connectivity/centres_orig.txt", usecols=[1, 2, 3], dtype=numpy.float64)
    centre_dict = {}
    for key, pos in zip(centre_keys, centre_positions):
        centre_dict[key] = pos

    result = []
    for idx, label in enumerate(centers_ordered):
        # if idx == 41:
        #     result.append(["unk_L", 0, 0, 0])
        position = centre_dict[label]
        result.append([label] + list(position))
    # result.append(["unk_R", 0, 0, 0])
    with open("Connectivity/centres.txt", "w") as outfile:
        outfile.write("\n".join(' '.join(str(v) for v in line) for line in result))


def create_zips():
    with ZipFile('connectivity_84.zip', 'w') as zip_object:
        zip_object.write('Connectivity/tract_lengths.txt')
        zip_object.write('Connectivity/weights.txt')
        zip_object.write('Connectivity/centres.txt')

    with ZipFile('surface_147k.zip', 'w') as zip_object:
        zip_object.write('Surface/vertices.txt')
        zip_object.write('Surface/triangles.txt')


if __name__ == "__main__":
    prepare_rm()
    prepare_centres()
    create_zips()
