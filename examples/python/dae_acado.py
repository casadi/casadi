#
#     This file is part of CasADi.
# 
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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
from numpy import *
import numpy as NP
import matplotlib.pyplot as plt

# CasADi
from casadi import *

try:
  ACADO_FCN_NUM_IN
except:
  print "This example needs ACADO"
  import sys
  sys.exit(0)
  
# Variables
t = ssym("t")
x = ssym("x")
l = ssym("l")
z = ssym("z")
u = ssym("u")

# Differential equation
f = vertcat([
  -x + 0.5*x*x + u + 0.5*z,
  x*x + 3.0*u*u,
  z + exp(z) - 1.0 + x
  ])

# The right hand side of the ACADO functions
acado_in = ACADO_FCN_NUM_IN * [[]]
acado_in[ACADO_FCN_T] = t     # Time
acado_in[ACADO_FCN_XD] = vertcat((x,l))  # Differential states
acado_in[ACADO_FCN_XA] = z    # Algebraic state
acado_in[ACADO_FCN_U] = u     # Control

# The DAE function
ffcn = SXFunction(acado_in,daeOut(f))

## Objective function
mfcn = SXFunction(acado_in,[l])

# Create ACADO solver
ocp_solver = AcadoOCP(ffcn,mfcn)

# Set options
t0 = 0.0
tf = 5.0
ocp_solver.setOption("start_time",t0)
ocp_solver.setOption("final_time",tf)

num_nodes = 20
ocp_solver.setOption("number_of_shooting_nodes",num_nodes)
#ocp_solver.setOption("max_num_iterations",30)
#ocp_solver.setOption("max_num_integrator_steps",10000)
#ocp_solver.setOption("dynamic_sensitivity","forward_sensitivities")
ocp_solver.setOption("kkt_tolerance",1e-5)
#ocp_solver.setOption("absolute_tolerance",1e-4)
#ocp_solver.setOption("integrator_tolerance",1e-8)
#ocp_solver.setOption("auto_init",True)
#ocp_solver.setOption("relaxation_parameter",1.5)

# Initialize
ocp_solver.init()

# Pass bounds
lbx0 = [1.0, 0.0, -NP.inf]
ubx0 = [1.0, 0.0, NP.inf]
ocp_solver.setInput(lbx0,"lbx0")
ocp_solver.setInput(ubx0,"ubx0")

# Solve
ocp_solver.solve()

# Time grid
t_opt = NP.linspace(t0,tf,num_nodes+1)

# Plot optimal state trajectory
x_opt = ocp_solver.output("x_opt")
x_opt = array(x_opt) # create numpy array
x_opt = x_opt.reshape(num_nodes+1, 3)

# Plot optimal control
u_opt = ocp_solver.output("u_opt")

# State derivatives
xdot = ssym("xdot")
ldot = ssym("ldot")

# The residual of the IDAS dae
dae_in = daeIn(t=t, x=vertcat((x,l)), z=z, xdot=vertcat((xdot,ldot)), p=u)

# The DAE residual
ode_res = vertcat((f[0]-xdot,f[1]-ldot))
alg_res = f[2]

# The DAE residual function
dae = SXFunction(dae_in,daeOut(ode=ode_res,alg=alg_res))

# Create an integrator
integrator = IdasIntegrator(dae)

# Set options
integrator.setOption("abstol",1e-6)
integrator.setOption("reltol",1e-6)
#integrator.setOption("exact_jacobian",True)
integrator.setOption("stop_at_end",True)
#integrator.setOption("suppress_algebraic",True)
#integrator.setOption("linear_solver_type","dense")
integrator.setOption("linear_solver_type","iterative")
integrator.setOption("t0",t_opt[0])
integrator.setOption("tf",t_opt[1])

# Initialize the integrator
integrator.init()

# Initial state
xi = (1.0,0.0)

# Simulated trajectory
x_opt2 = array((num_nodes+1)*[xi])

# Loop over the time points
for i in range(num_nodes):
  
  # Set the control
  ui = u_opt[i]
  integrator.setInput(ui,"p")

  # Pass the state
  integrator.setInput(xi,"x0")
  
  # Integrate
  integrator.evaluate()
  
  # Get the state
  xi = integrator.output("xf")
  x_opt2[i+1] = xi.data()

plt.figure(1)
plt.clf()
plt.hold(True)
plt.plot(t_opt,x_opt[:,0],'*')
plt.plot(t_opt,x_opt2[:,0],'-')
plt.title("DIFFERENTIAL STATE  x")
plt.legend(('ACADO','IDAS'))

plt.figure(2)
plt.clf()
plt.hold(True)
plt.plot(t_opt,x_opt[:,2],'*')
#plt.plot(t_opt,x_opt2[:,2],'-')
plt.title("ALGEBRAIC STATE  z")
#plt.legend(('ACADO','IDAS'))

plt.figure(3)
plt.clf()
plt.hold(True)
plt.plot(t_opt,x_opt[:,1],'*')
plt.plot(t_opt,x_opt2[:,1],'-')
plt.title("Lagrange term l")
plt.legend(('ACADO','IDAS'))

plt.figure(3)
plt.clf()
plt.hold(True)
plt.plot(t_opt,x_opt[:,1],'*')
plt.plot(t_opt,x_opt2[:,1],'-')
plt.title("Lagrange term l")
plt.legend(('ACADO','IDAS'))

plt.figure(4)
plt.clf()
plt.hold(True)
plt.plot(t_opt,u_opt[:],'--')
#plt.plot(t_opt,u_opt2[:],'-')
plt.title("CONTROL u")
#plt.legend(('ACADO','IDAS'))

plt.show()
