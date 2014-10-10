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
from pymodelica import compile_jmu
from pyjmi import JMUModel
import pymodelica
import os, zipfile

# Compile the modelica source
curr_dir = os.path.dirname(os.path.abspath(__file__))
jmu_name = compile_jmu("TimedVariablesTest", curr_dir+"/path_constraints.mop",'optimica','ipopt',{'generate_xml_equations':True, 'generate_fmi_me_xml':False})
sfile = zipfile.ZipFile(curr_dir+'/TimedVariablesTest.jmu','r')
mfile = sfile.extract('modelDescription.xml','.')
os.remove(curr_dir+'/TimedVariablesTest.jmu')

# Import into CasADi
ocp = casadi.SymbolicOCP()
ocp.parseFMI('modelDescription.xml',{'sort_equations':False,'eliminate_dependent':False})
ocp.sortType(True) # temporary solution: enables the new sorting
print ocp

