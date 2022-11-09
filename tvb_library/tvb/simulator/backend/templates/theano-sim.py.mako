# -*- coding: utf-8 -*-
#
#
# TheVirtualBrain-Scientific Package. This package holds all simulators, and
# analysers necessary to run brain-simulations. You can use it stand alone or
# in conjunction with TheVirtualBrain-Framework Package. See content of the
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

import numpy as np
import theano
import theano.tensor as tt

<%include file="theano-coupling.py.mako" />
<%include file="theano-dfuns.py.mako" />
<%include file="theano-integrate.py.mako" />

<%
    from tvb.simulator.integrators import IntegratorStochastic
    stochastic = isinstance(sim.integrator, IntegratorStochastic)
    any_delays = sim.connectivity.idelays.any()
%>

def kernel(state, weights, trace, parmat
           ${', nsig' if stochastic else ''}
           ${', idelays' if any_delays else ''}
           ):

    # problem dimensions
    n_node = ${sim.connectivity.weights.shape[0]}
    n_svar = ${len(sim.model.state_variables)}
    n_cvar = ${len(sim.model.cvar)}
    nt = ${int(sim.simulation_length/sim.integrator.dt)}

    # work space arrays
    dX = tt.zeros((${sim.integrator.n_dx}, n_svar, n_node))
    cX = tt.zeros((n_cvar, n_node))

    for t in range(nt):
        state = integrate(state, weights, parmat, dX, cX
           ${', nsig' if stochastic else ''}
           ${', idelays' if any_delays else ''}
           )
        trace = tt.set_subtensor(trace[t], state[:,0])

    return trace
