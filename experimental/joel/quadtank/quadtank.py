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
import scipy.integrate as integr

# JModelica
from jmodelica.jmi import compile_jmu
from jmodelica.jmi import JMUModel
import jmodelica

# CasADi
from casadi import *

curr_dir = os.path.dirname(os.path.abspath(__file__))

try:
  # Try the old Jmodelica syntax
  jmu_name = compile_jmu("QuadTank_pack.QuadTank_Opt", curr_dir+"/QuadTank.mop",'optimica','ipopt',{'generate_xml_equations':True, 'generate_fmi_xml':True})
except jmodelica.compiler.UnknownOptionError:
  # Try the new jmodelica syntax
  jmu_name = compile_jmu("QuadTank_pack.QuadTank_Opt", curr_dir+"/QuadTank.mop",'optimica','ipopt',{'generate_xml_equations':True, 'generate_fmi_me_xml':True})
  
import zipfile
sfile = zipfile.ZipFile(curr_dir+'/QuadTank_pack_QuadTank_Opt.jmu','r')
mfile = sfile.extract('modelDescription.xml','.')

# Allocate a parser and load the xml
parser = FMIParser('modelDescription.xml')

# Dump representation to screen
print "XML representation"
print parser

# Obtain the symbolic representation of the OCP
ocp = parser.parse()

# Sort the variables according to BLT
ocp.sortBLT()

# Print the ocp to screen
print ocp

# Make explicit
ocp.makeExplicit()
  
# Print again
print ocp

# The right hand side of the ACADO functions
acado_in = ACADO_FCN_NUM_IN * [[]]

# Time
acado_in[ACADO_FCN_T] = ocp.t

# Differential state
acado_in[ACADO_FCN_XD] = ocp.xd

# Algebraic state
acado_in[ACADO_FCN_XA] = ocp.xa

# Control
acado_in[ACADO_FCN_U] = ocp.u

# Parameter
acado_in[ACADO_FCN_P] = ocp.p

# The DAE function
ffcn_out = list(ocp.diffeq) + list(ocp.algeq)
ffcn = SXFunction(acado_in,[ffcn_out])

ffcn.init()

def res2(y,t):
  ffcn.setInput(y, "fcn_xd");
  ffcn.evaluate()
  return list(ffcn.getOutputData())

x_0 = [0.0, 0.01, 0.01, 0.01, 0.01]
t_sim = NP.linspace(0.,2000.,500)

u_A = [2.,2]
ffcn.setInput(u_A, "fcn_u");
y_sim = integr.odeint(res2,x_0,t_sim)

# Plot
plt.figure(1)
plt.clf()
plt.subplot(211)
plt.plot(t_sim,y_sim[:,1:5])
plt.grid()

u_B = [2.5,2.5]
ffcn.setInput(u_B, "fcn_u");
y_sim = integr.odeint(res2,x_0,t_sim)

# Plot
plt.figure(1)
plt.subplot(212)
plt.plot(t_sim,y_sim[:,1:5])
plt.grid()

# Objective function
mfcn = SXFunction(acado_in,[ocp.mterm])

# Initial constraint function
rfcn = SXFunction(acado_in,[ocp.initeq])

# Create ACADO solver
ocp_solver = AcadoInterface(ffcn,mfcn,Function(),rfcn)

# Set options
ocp_solver.setOption("start_time",ocp.t0)
ocp_solver.setOption("final_time",ocp.tf)
num_nodes = 20
ocp_solver.setOption("number_of_shooting_nodes",num_nodes)
ocp_solver.setOption("max_num_iterations",200)
ocp_solver.setOption("kkt_tolerance",1e-6)

# Initialize
ocp_solver.init()

# Set bounds on states
x_guess = (num_nodes+1) * [0, 0.1, 0.1, 0.1, 0.1]
u0 = (num_nodes+1) * u_A
ocp_solver.setInput(x_guess,"x_guess")
ocp_solver.setInput(u0,"u_guess")
  
ocp_solver.setInput([2,2],"lbu")
ocp_solver.setInput([10,10],"ubu")
ocp_solver.setInput((num_nodes+1)*[5,5],"u_guess")

# Solve the optimal control problem
ocp_solver.solve()

# Optimal cost
cost = ocp_solver.getOutputData(ACADO_COST)[0]
print "optimal cost = ", cost

# Print optimal parameters
popt = ocp_solver.getOutputData(ACADO_P_OPT)
print "optimal parameter values = ", popt

# Time grid
t_opt = linspace(0,ocp.tf,num_nodes+1)

# Plot optimal control
u_opt = ocp_solver.getOutputData(ACADO_U_OPT)
u1 = []
for i in range(0,len(u_opt),2):
  u1.append(u_opt[i])

u2 = []
for i in range(1,len(u_opt),2):
  u2.append(u_opt[i])
  
plt.figure(2)
plt.subplot(211)
plt.plot(t_opt,u1)
plt.subplot(212)
plt.plot(t_opt,u2)

# Plot optimal state trajectory
x_opt = ocp_solver.getOutputData(ACADO_X_OPT)
x_opt = array(x_opt) # create numpy array
x_opt = x_opt.reshape(num_nodes+1, 5)
plt.figure(3)
plt.plot(t_opt,x_opt[:,0])

plt.figure(4)
plt.plot(t_opt,x_opt[:,1:5])

plt.show()

