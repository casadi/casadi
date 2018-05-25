#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             K.U. Leuven. All rights reserved.
#     Copyright (C) 2011-2014 Greg Horn
#
#     CasADi is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 3 of the License, or (at your option) any later version.
#
#     CasADi is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public
#     License along with CasADi; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
#
from casadi import *

"""
This example demonstrates how NL-files, which can be generated
by AMPl or Pyomo, can be imported in CasADi and solved using
e.g. the interface to AMPL

Joel Andersson
2012

"""

# Create an NLP instance
nl = NlpBuilder()

# Parse an NL-file
nl.import_nl("../nl_files/hs107.nl",{"verbose":False})

# NLP solver options
opts = {}
opts["expand"] = True
# opts["verbose"] = True
# opts["ipopt"] = dict(max_iter=10, linear_solver="ma57", hessian_approximation="limited-memory")

# Create an NLP solver
nlpsol = nlpsol("nlpsol", "ipopt", nl, opts)

# Solve NLP
res = nlpsol(lbx=nl.x_lb,
             ubx=nl.x_ub,
             lbg=nl.g_lb,
             ubg=nl.g_ub,
             x0=nl.x_init)
