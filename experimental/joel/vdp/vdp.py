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
# -*- coding: utf-8 -*-
import os
import sys
from numpy import *
import numpy as NP
import matplotlib.pyplot as plt
import zipfile

# JModelica
from jmodelica.jmi import compile_jmu
from jmodelica.jmi import JMUModel
import jmodelica

# CasADi
from casadi import *

curr_dir = os.path.dirname(os.path.abspath(__file__));
try:
  # Try the old Jmodelica syntax
  jmu_name = compile_jmu("VDP_pack.VDP_Opt", curr_dir+"/VDP.mop",'optimica','ipopt',{'generate_xml_equations':True, 'generate_fmi_xml':False})
except jmodelica.compiler.UnknownOptionError:
  # Try the new jmodelica syntax
  jmu_name = compile_jmu("VDP_pack.VDP_Opt", curr_dir+"/VDP.mop",'optimica','ipopt',{'generate_xml_equations':True, 'generate_fmi_me_xml':False})
  
if True:
  vdp = JMUModel(jmu_name)
  res = vdp.optimize()

  # Extract variable profiles
  x1=res['x1']
  x2=res['x2']
  u=res['u']
  t=res['time']
  cost=res['cost']

  # Plot
  plt.figure(1)
  plt.clf()
  plt.subplot(311)
  plt.plot(t,x1)
  plt.grid()
  plt.ylabel('x1')
        
  plt.subplot(312)
  plt.plot(t,x2)
  plt.grid()
  plt.ylabel('x2')
        
  plt.subplot(313)
  plt.plot(t,u)
  plt.grid()
  plt.ylabel('u')
  plt.xlabel('time')

sfile = zipfile.ZipFile(curr_dir+'/VDP_pack_VDP_Opt.jmu','r')
mfile = sfile.extract('modelDescription.xml','.')
os.remove('VDP_pack_VDP_Opt.jmu')
os.rename('modelDescription.xml','vdp.xml')

# Allocate a parser and load the xml
parser = FMIParser('vdp.xml')

# Dump representation to screen
print "XML representation"
print parser

# Obtain the symbolic representation of the OCP
ocp = parser.parse()

# Print the ocp to screen
print ocp

# Sort the variables according to type
var = OCPVariables(ocp.variables)

# The right hand side of the ACADO functions
acado_in = ACADO_FCN_NUM_IN * [[]]

# Time
acado_in[ACADO_FCN_T] = [var.t_]

# Convert stl vector of variables to list of expressions
def toList(v, der=False):
  ret = []
  for i in v:
    if der:
      ret.append(i.der())
    else:
      ret.append(i.var())
  return ret

# Differential state
acado_in[ACADO_FCN_XD]  = toList(ocp.x_)

# Algebraic state
acado_in[ACADO_FCN_XA] = toList(ocp.z_)

# Control
acado_in[ACADO_FCN_U] = toList(ocp.u_)

# Parameter
acado_in[ACADO_FCN_P] = toList(ocp.p_)

# State derivative
acado_in[ACADO_FCN_XDOT] = toList(ocp.x_,True)

# The DAE function
ffcn_out = list(ocp.dae) + list(ocp.ae)

ffcn = SXFunction(acado_in,[ffcn_out])

# Objective function
mfcn = SXFunction(acado_in,[ocp.mterm])

# Path constraint function
cfcn = SXFunction(acado_in,[ocp.cfcn])
  
# Initial constraint function
rfcn = SXFunction(acado_in,[ocp.initeq])

# Create ACADO solver
ocp_solver = AcadoInterface(ffcn,mfcn,cfcn,rfcn)

# Create an integrator
dae_in = DAE_NUM_IN * [[]]
dae_in[DAE_T] = acado_in[ACADO_FCN_T]
dae_in[DAE_Y] = acado_in[ACADO_FCN_XD] + acado_in[ACADO_FCN_XA]
dae_in[DAE_YDOT] = acado_in[ACADO_FCN_XDOT] + list(ssym("zdot",len(acado_in[ACADO_FCN_XA])))
dae_in[DAE_P] = acado_in[ACADO_FCN_P] + acado_in[ACADO_FCN_U]
dae = SXFunction(dae_in,[ffcn_out])

integrator = IdasIntegrator(dae)
#integrator.setOption("exact_jacobian",True)
#integrator.setOption("linear_multistep_method","bdf") # adams or bdf
#integrator.setOption("nonlinear_solver_iteration","newton") # newton or functional
integrator.setOption("number_of_fwd_dir",4)
integrator.setOption("number_of_adj_dir",0)
integrator.setOption("fsens_err_con",True)
integrator.setOption("quad_err_con",True)
integrator.setOption("abstol",1e-8)
integrator.setOption("reltol",1e-8)
integrator.setOption("is_differential",len(acado_in[ACADO_FCN_XD])*[1] + len(acado_in[ACADO_FCN_XA])*[0])


# Pass the integrator to ACADO
ocp_solver.setIntegrator(integrator)

# Set options
ocp_solver.setOption("start_time",ocp.t0)
ocp_solver.setOption("final_time",ocp.tf)
num_nodes = 30
ocp_solver.setOption("number_of_shooting_nodes",num_nodes)
ocp_solver.setOption("max_num_iterations",100)
ocp_solver.setOption("kkt_tolerance",1e-4)
ocp_solver.setOption("integrator","casadi")
ocp_solver.setOption("integrator_tolerance",1e-6)

# Initialize
ocp_solver.init()

# Set bounds on states
cfcn_lb = []
for i in ocp.cfcn_lb:
  cfcn_lb.append(float(i))
ocp_solver.setInput(cfcn_lb,"lbc")
  
cfcn_ub = []
for i in ocp.cfcn_ub:
  cfcn_ub.append(float(i))
ocp_solver.setInput(cfcn_ub,"ubc")
  
# Solve the optimal control problem
ocp_solver.solve()

# Print optimal cost
cost = ocp_solver.getOutputData(ACADO_COST)[0]
print "optimal cost = ", cost

# Print optimal parameters
popt = ocp_solver.getOutputData(ACADO_P_OPT)
print "optimal parameter values = ", popt

# Time grid
t_opt = NP.linspace(0,ocp.tf,num_nodes+1)

# Plot optimal control
u_opt = ocp_solver.getOutputData(ACADO_U_OPT)
plt.figure(3)
plt.plot(t_opt,u_opt)

# Plot optimal state trajectory
x_opt = ocp_solver.getOutputData(ACADO_X_OPT)
x_opt = array(x_opt) # create numpy array
x_opt = x_opt.reshape(num_nodes+1, 3)
plt.figure(4)
plt.plot(t_opt,x_opt)

# Show the plots
plt.ion()
plt.show()
  
