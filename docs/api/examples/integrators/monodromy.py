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
#! Monodromy matrix
#! =====================
from casadi import *
from numpy import *
from pylab import *

#! We will investigate the monodromy matrix with the help of a simple 2-state system, as found in 1. Nayfeh AH, Balachandran B. Applied nonlinear dynamics. 1995. Available at: http://onlinelibrary.wiley.com/doi/10.1002/9783527617548.biblio/summary [Accessed June 16, 2011], page 52.
#$ $\dot{x_1} = x_2$
#$ $\dot{x_2} = -(-w_0^2 x_1 + a_3 x_1^3 + a_5 x_1^5) - (2 mu_1 x_2 + mu_3 x_2^3) + f$.

x  = SX.sym("x",2)
x1 = x[0]
x2 = x[1]

w0 = SX.sym("w0")
a3 = SX.sym("a3")
a5 = SX.sym("a5")
mu1 = SX.sym("mu1")
mu3 = SX.sym("mu3")
ff = SX.sym("f")

tf = 40

params = vertcat([w0,a3,a5,mu1,mu3,ff])
rhs    = vertcat([x2,(-(-w0**2 *x1 + a3*x1**3 + a5*x1**5) - (2 *mu1 *x2 + mu3 * x2**3))/100+ff])

f=SXFunction(daeIn(x=x,p=params),daeOut(ode=rhs))
f.init()

t = SX.sym("t")
cf=SXFunction(controldaeIn(t=t, x=x, p=vertcat([w0,a3,a5,mu1,mu3]), u=ff),[rhs])
cf.init()

integrator = Integrator("cvodes", f)
integrator.setOption("tf",tf)
integrator.setOption("reltol",1e-10)
integrator.setOption("abstol",1e-10)
integrator.setOption("fsens_err_con",True)
integrator.init()

N = 500

#! Let's get acquainted with the system by drawing a phase portrait
ts = linspace(0,tf,N)

sim = Simulator(integrator,ts)
sim.init()

w0_ = 5.278
params_ = [ w0_, -1.402*w0_**2,  0.271*w0_**2,0,0,0 ]

sim.setInput(params_,"p")

x2_0 = 0
figure(1)
for x1_0 in [-3.5,-3.1,-3,-2,-1,0]:
  sim.setInput([x1_0,x2_0],"x0")
  sim.evaluate()
  plot(sim.getOutput()[0,:],sim.getOutput()[1,:],'k')

title('phase portrait for mu_1 = 0, mu_2 = 0')
xlabel('x_1')
ylabel('x_2')

show()

x0 = DMatrix([-3.1,0])

#! Monodromy matrix at tf - Jacobian of integrator
#! ===============================================
#! First argument is input index, second argument is output index
jac = integrator.jacobian("x0","xf")
jac.init()

jac.setInput(x0,"x0")
jac.setInput(params_,"p")
jac.evaluate()

Ji = jac.getOutput()

print Ji

#! Note the remainder of this file depends on Jacobian of Simulator,
#!  a feature that is currently disabled 
#! https://github.com/casadi/casadi/issues/929

# #! Monodromy matrix at various instances - Jacobian of Simulator
# #! =============================================================

# jacsim = sim.jacobian(INTEGRATOR_X0,0)
# jacsim.init()

# jacsim.setInput(x0,"x0")
# jacsim.setInput(params_,"p")
# jacsim.evaluate()

# #! For each of the 500 intervals, we have a 2-by-2 matrix as output
# print "jacsim.output().shape = ", jacsim.output().shape

# #! Show only the last 3 intervals.
# print jacsim.getOutput()[-3*2:,:]

# Js = jacsim.getOutput()[-2:,:]


# e = max(fabs(Js - Ji))/max(fabs(Js))

# # Assert that the two methods yield identical results
# assert(e < 1e-6)

# #! Monodromy matrix at various instances - Jacobian of ControlSimulator
# #! ====================================================================

# csim = ControlSimulator(cf,linspace(0,tf,50))
# csim.setOption("nf",10)
# csim.setOption("integrator","cvodes")
# csim.setOption("integrator_options",{"reltol":1e-11,"abstol":1e-11, "fsens_err_con": True})
# csim.init()

# jaccsim = csim.jacobian(CONTROLSIMULATOR_X0,0)
# jaccsim.init()
# jaccsim.setInput(params_[:-1],"p")
# jaccsim.setInput(x0,"x0")
# jaccsim.setInput(0,"u")
# jaccsim.evaluate()

# #! For each of the 500 intervals, we have a 2-by-2 matrix as output
# print "jaccsim.output().shape = ", jaccsim.output().shape

# #! Show only the last 3 intervals.
# print jaccsim.getOutput()[-3*2:,:]
# Jcs = jaccsim.getOutput()[-2:,:]

# e = max(fabs(Jcs - Js))/max(fabs(Js))

# # Assert that the two methods yield identical results
# assert(e < 1e-5)

# #! Intuitive interpretation
# #! ========================

# sim.setInput(x0,"x0")
# sim.setInput(params_,"p")
# sim.evaluate()
# unperturbed_output = sim.getOutput()

# circle = array([[sin(x),cos(x)] for x in numpy.linspace(-pi/2,3/2.0*pi,100)]).T
# circle = hstack((circle,circle[:,50:51]))



# for t in range(0,N/5,2):
#   J = jacsim.getOutput()[t*2:(t+1)*2,:]
#   if t < 10:
#     scale = 0.1
#   else:
#     scale = 0.01
#   e=scale*mul(J,circle).T
#   e[:,0] += sim.getOutput()[t,0]
#   e[:,1] += sim.getOutput()[t,1]
#   if t < 10 :
#     plot(e[:,0],e[:,1],color='red')
#   else:
#     plot(e[:,0],e[:,1],color='blue')
    
# show()
# #! Consider the case of perturbation simulation with a slightly perturbed initial condition

# sim.setInput(x0,"x0")
# sim.setInput(params_,"p")
# sim.evaluate()
# unperturbed_output = sim.getOutput()

# perturb = DMatrix([1e-2,0])
# sim.setInput(x0+perturb,"x0")
# sim.setInput(params_,"p")
# sim.evaluate()
# perturbed_output = sim.getOutput()

# figure(2)

# title('Evolution of a perturbation')
# plot(ts,perturbed_output-unperturbed_output)

# effects = DMatrix.zeros(N,2)

# for t in range(N):
#   effects[t,:] = mul(jacsim.getOutput()[t*2:(t+1)*2,:],perturb).T
  
# plot(ts,effects)

# legend(('x_1','x_2','perturbed(x_1)','preturbed(y_2)'))
# xlabel('t')

# show()

# figure(3)
# linear_perturbed = unperturbed_output.reshape((2*N,1)) + mul(jacsim.getOutput(),perturb)

# title('phase portrait perturbation')
# plot(unperturbed_output[:,0],unperturbed_output[:,1])
# plot(perturbed_output[:,0],perturbed_output[:,1])
# plot(linear_perturbed[0:N/2:2],linear_perturbed[1:N/2:2])

# legend(('nominal','pertubed','monodromy prediction'))


show()

