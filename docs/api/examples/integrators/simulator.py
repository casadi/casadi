#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             KU Leuven. All rights reserved.
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
# Simulator
# =====================

import casadi as ca
from numpy import *
from pylab import *

# We will investigate the working of Simulator with the help of the parametrically exited Duffing equation:
#
# $\ddot{u}+\dot{u}-\epsilon (2 \mu \dot{u}+\alpha u^3+2 k u \cos(\Omega t))$ with $\Omega = 2 + \epsilon \sigma$.

t = ca.SX.sym('t')

u = ca.SX.sym('u') 
v = ca.SX.sym('v') 
states = ca.vertcat(u,v)

eps   = ca.SX.sym('eps')
mu    = ca.SX.sym('mu')
alpha = ca.SX.sym('alpha')
k     = ca.SX.sym('k')
sigma = ca.SX.sym('sigma')
Omega = 2 + eps*sigma

params = ca.vertcat(eps,mu,alpha,k,sigma)
rhs    = ca.vertcat(v,-u-eps*(2*mu*v+alpha*u**3+2*k*u*cos(Omega*t)))

# We will simulate over 50 seconds, 1000 timesteps.

dae={'x':states, 'p':params, 't':t, 'ode':rhs}
ts = linspace(0, 50, 1000)
integrator = ca.integrator('integrator', 'cvodes', dae, 0, ts)

sol = integrator(x0=[1,0], p=[0.1,0.1,0.1,0.3,0.1])

# Plot the solution

plot(array(sol['xf'])[0,:], array(sol['xf'])[1,:])
xlabel('u')
ylabel('u_dot')
show()

