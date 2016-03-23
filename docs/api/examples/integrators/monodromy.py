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

params = vertcat(w0,a3,a5,mu1,mu3,ff)
rhs    = vertcat(x2,(-(-w0**2 *x1 + a3*x1**3 + a5*x1**5) - (2 *mu1 *x2 + mu3 * x2**3))/100+ff)

dae={'x':x, 'p':params, 'ode':rhs}

# t = SX.sym("t")
# cf=SX.fun("cf", controldaeIn(t=t, x=x, p=vertcat(w0,a3,a5,mu1,mu3), u=ff),[rhs])

# opts = {}
# opts["tf"] = tf
# opts["reltol"] = 1e-10
# opts["abstol"] = 1e-10
# opts["fsens_err_con"] = True
# integrator = integrator("integrator", "cvodes", dae, opts)

# N = 500

# #! Let's get acquainted with the system by drawing a phase portrait
# ts = linspace(0,tf,N)

# sim = Simulator("sim", integrator,ts)

# w0_ = 5.278
# params_ = [ w0_, -1.402*w0_**2,  0.271*w0_**2,0,0,0 ]

# sim.setInput(params_,"p")

# x2_0 = 0
# figure(1)
# for x1_0 in [-3.5,-3.1,-3,-2,-1,0]:
#   sim.setInput([x1_0,x2_0],"x0")
#   sim.evaluate()
#   plot(sim.getOutput()[0,:],sim.getOutput()[1,:],'k')

# title('phase portrait for mu_1 = 0, mu_2 = 0')
# xlabel('x_1')
# ylabel('x_2')

# show()

# x0 = DM([-3.1,0])

# #! Monodromy matrix at tf - Jacobian of integrator
# #! ===============================================
# #! First argument is input index, second argument is output index
# jac = integrator.jacobian("x0","xf")

# jac.setInput(x0,"x0")
# jac.setInput(params_,"p")
# jac.evaluate()

# Ji = jac.getOutput()

# print Ji
