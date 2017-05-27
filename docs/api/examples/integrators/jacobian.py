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
#! Integrator jacobian
#! =====================
from casadi import *
from numpy import *

#! We will investigate the working of integrator jacobian with the help of the parametrically exited Duffing equation:
#!
#$ $\ddot{u}+\dot{u}-\epsilon (2 \mu \dot{u}+\alpha u^3+2 k u \cos(\Omega t))$ with $\Omega = 2 + \epsilon \sigma$.

t = SX.sym('t')
u = SX.sym('u')
v = SX.sym('v')

eps   = SX.sym('eps')
mu    = SX.sym('mu')
alpha = SX.sym('alpha')
k     = SX.sym('k')
sigma = SX.sym('sigma')
Omega = 2 + eps*sigma

params = vertcat(eps,mu,alpha,k,sigma)
states = vertcat(u,v)
rhs    = vertcat(v,-u-eps*(2*mu*v+alpha*u**3+2*k*u*cos(Omega*t)))

dae = {'x':states, 'p':params, 't':t, 'ode':rhs}

F = integrator('F', 'cvodes', dae)

#! First argument is input index, secpnd argument is output index
jac = F.factory('jac_F', F.name_in(), ['jac:xf:p'])
