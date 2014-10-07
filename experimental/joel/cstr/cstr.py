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
import numpy as NP
from numpy import *
import matplotlib.pyplot as plt
import zipfile
import time

# JModelica
from jmodelica.jmi import compile_jmu
from jmodelica.jmi import JMUModel
import jmodelica

# CasADi
from casadi import *

curr_dir = os.path.dirname(os.path.abspath(__file__))

try:
  # Try the old Jmodelica syntax
  jmu_name = compile_jmu("CSTR.CSTR_Opt", curr_dir+"/CSTR.mop",'optimica','ipopt',{'generate_xml_equations':True, 'generate_fmi_xml':False})
except jmodelica.compiler.UnknownOptionError:
  # Try the new jmodelica syntax
  jmu_name = compile_jmu("CSTR.CSTR_Opt", curr_dir+"/CSTR.mop",'optimica','ipopt',{'generate_xml_equations':True, 'generate_fmi_me_xml':False})

Tc_ref = 280
c_ref = 338.775766
T_ref = 280.099198

c_init = 956.271065
T_init = 250.051971

sfile = zipfile.ZipFile(curr_dir+'/CSTR_CSTR_Opt.jmu','r')
mfile = sfile.extract('modelDescription.xml','.')

# Allocate a parser and load the xml
parser = FMIParser('modelDescription.xml')

# Dump representation to screen
print "XML representation"
print parser

# Obtain the symbolic representation of the OCP
ocp = parser.parse()

# Print the ocp to screen
print ocp

# Convert stl vector of variables to an array of expressions
def toArray(v, der=False):
  ret = []
  for i in v:
    if der:
      ret.append(i.der())
    else:
      ret.append(i.var())
  return array(ret,dtype=SX)

# Variables
t = ocp.t_
x = toArray(ocp.xd_)
xdot = toArray(ocp.xd_,True)
xa = toArray(ocp.xa_)
p = toArray(ocp.p_)
u = toArray(ocp.u_)

# Get scaling factors values
def getNominal(v):
  ret = []
  for i in v:
    ret.append(i.getNominal())
  return array(ret,dtype=SX)

x_sca = getNominal(ocp.xd_)
xa_sca = getNominal(ocp.xa_)
u_sca = getNominal(ocp.u_)
p_sca = getNominal(ocp.p_)

print "x_sca = ", repr(x_sca)
print "xa_sca = ", repr(xa_sca)
print "u_sca = ", repr(u_sca)
print "p_sca = ", repr(p_sca)

# The old variables expressed in the normalized variables
x_old = x_sca*x
xdot_old = x_sca*xdot
xa_old = xa_sca*xa
u_old = u_sca*u
p_old = p_sca*p

# scale a function
def scale_exp(f_old):
  ffcn_old = SXFunction([x,xdot,xa,p,u],[f_old])
  ffcn_old.init()
  f_new = ffcn_old.eval([x_old,xdot_old,xa_old,p_old,u_old])[0]
  return array(f_new) # NB! typemap will change

# Create an integrator
dae_in = DAE_NUM_IN * [[]]
dae_in[DAE_T] = [t]
dae_in[DAE_Y] = x
dae_in[DAE_YDOT] = xdot
dae_in[DAE_Z] = xa
dae_in[DAE_P] = concatenate((p,u))
dae = SXFunction(dae_in,[scale_exp(ocp.dynamic_eq_)/[[1e7], [1000.], [350.]]])

# Number of shooting nodes
num_nodes = 100

# Bounds on states
cfcn_lb = []
for i in ocp.cfcn_lb:
  cfcn_lb.append(float(i))  
cfcn_ub = []
for i in ocp.cfcn_ub:
  cfcn_ub.append(float(i))

# Initial guess for the state
x0 = array([0.0, T_init, c_init])/x_sca

# Initial guess for the control
u0 = array([280])/u_sca

# Create integrators
integrators = []
for i in range(num_nodes*2):
  integrator = IdasIntegrator(dae)
  integrator.setOption("number_of_fwd_dir",1)
  integrator.setOption("number_of_adj_dir",0)
  integrator.setOption("fsens_err_con",True)
  integrator.setOption("quad_err_con",True)
  integrator.setOption("abstol",1e-8)
  integrator.setOption("reltol",1e-8)
  integrator.init()
  integrators.append(integrator)

# Number of differential states
nx = 3

# Number of controls
nu = 1

# Create a multiple shooting discretization
#ms = MultipleShooting(integrator,num_nodes,nx,nu)

## Copy data
#ms.tf_ = ocp.tf

#ms.u_init_  = DVector([280]/u_sca)
#ms.u_min_   = DVector([230]/u_sca)
#ms.u_max_   = DVector([370]/u_sca)

#ms.x_init_  = DVector(x0)
#ms.x_min_   = DVector([-inf,-inf,-inf]/x_sca)
#ms.x_max_   = DVector([inf, 350, inf]/x_sca)
#ms.xf_min_  = DVector([-inf,-inf,-inf]/x_sca)
#ms.xf_max_  = DVector([inf, 350,  inf]/x_sca)
#ms.x0_min_  = DVector([0,250.052,956.271]/x_sca)
#ms.x0_max_  = DVector([0,250.052,956.271]/x_sca)

#ms.init()
  
#solver = IpoptSolver(ms.F_,ms.G_,Function(),ms.J_)
#solver.setOption("tol",1e-5)
#solver.setOption("hessian_approximation", "limited-memory")
#solver.setOption("max_iter",100)
#solver.setOption("linear_solver","ma57")
##  solver.setOption("derivative_test","first-order")

#solver.setOption("verbose",True)
#solver.init()

## Set bounds and initial guess
#solver.setInput(ms.V_min_,  "lbx")
#solver.setInput(ms.V_max_,  "ubx")
#solver.setInput(ms.V_init_,  "x0")

#solver.setInput(ms.G_min_,"lbg")
#solver.setInput(ms.G_max_,"ubg")

# Solve the problem
#solver.solve()



# The right hand side of the ACADO functions
acado_in = ACADO_FCN_NUM_IN * [[]]
acado_in[ACADO_FCN_T] = [t]
acado_in[ACADO_FCN_XD]  = x
acado_in[ACADO_FCN_XA] = z
acado_in[ACADO_FCN_U] = u
acado_in[ACADO_FCN_P] = p
acado_in[ACADO_FCN_XDOT] = xdot

# The DAE function
ffcn = SXFunction(acado_in,[scale_exp(ocp.dynamic_eq_)/[[1e7], [1000.], [350.]]])

# Objective function
mfcn = SXFunction(acado_in,[scale_exp(ocp.mterm)/1e7])

# Path constraint function
cfcn = SXFunction(acado_in,[scale_exp(ocp.cfcn)])
  
# Initial constraint function
rfcn = SXFunction(acado_in,[scale_exp(ocp.initial_eq_)])

# Create ACADO solver
ocp_solver = AcadoInterface(ffcn,mfcn,cfcn,rfcn)

# Pass the integrators to ACADO
ocp_solver.setIntegrators(integrators)

# Set options
ocp_solver.setOption("start_time",ocp.t0)
ocp_solver.setOption("final_time",ocp.tf)
ocp_solver.setOption("number_of_shooting_nodes",num_nodes)
ocp_solver.setOption("max_num_iterations",100)
ocp_solver.setOption("kkt_tolerance",6e-7)
#ocp_solver.setOption("integrator","casadi")
ocp_solver.setOption("integrator_tolerance",1e-8)

# Initialize
ocp_solver.init()

# Pass bounds and initial guess
ocp_solver.setInput(cfcn_lb,"lbc")
ocp_solver.setInput(cfcn_ub,"ubc")
ocp_solver.setInput((num_nodes+1)*list(x0),"x_guess")
ocp_solver.setInput((num_nodes+1)*list(u0),"u_guess")

# Solve the optimal control problem
t1 = time.time()
ocp_solver.solve()
t2 = time.time()

print 'Optimization took %0.3f ms' % ((t2-t1)*1000.0)

# Print optimal cost
cost = ocp_solver.getOutput("cost")[0]
print "optimal cost = ", cost

# Print optimal parameters
popt = ocp_solver.getOutput("p_opt")
print "optimal parameter values = ", popt

# Time grid
t_opt = NP.linspace(0,ocp.tf,num_nodes+1)

# Plot optimal control
u_opt = ocp_solver.getOutput("u_opt")
plt.figure(1)
plt.plot(t_opt,trans(u_opt))

# Plot optimal state trajectory
x_opt = ocp_solver.getOutput("x_opt")
x_opt = array(x_opt) # create numpy array
x_opt = x_opt.reshape(num_nodes+1, 3)
plt.figure(2)
plt.clf()

plt.subplot(3,1,1)
plt.plot(t_opt,x_opt[:,0])

plt.subplot(3,1,2)
plt.plot(t_opt,x_opt[:,1])

plt.subplot(3,1,3)
plt.plot(t_opt,x_opt[:,2])

print repr(t_opt)
print repr(x_opt)
print repr(u_opt)



# Show the plots
plt.ion()
plt.show()

