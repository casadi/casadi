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
#! NLPImplicitSolver
#! =====================
from casadi import *
from numpy import *
from pylab import *

#! We will investigate the working of rootfinder with the help of the parametrically exited Duffing equation.
#!
#$ $\ddot{u}+\dot{u}-\epsilon (2 \mu \dot{u}+\alpha u^3+2 k u \cos(\Omega t))$ with $\Omega = 2 + \epsilon \sigma$. \\
#$
#$ The first order solution is $u(t)=a \cos(\frac{1}{2} \Omega t-\frac{1}{2} \gamma)$ with the modulation equations: \\
#$ $\frac{da}{d\epsilon t} = - \left[ \mu a + \frac{1}{2} k a \sin \gamma \right]$ \\
#$ $a \frac{d\gamma}{d\epsilon t} = - \left[ -\sigma a + \frac{3}{4} \alpha a^3 + k a \cos \gamma \right]$ \\
#$
#$ We seek the stationair solution to these modulation equations.

#! Parameters
eps   = SX.sym("eps")
mu    = SX.sym("mu")
alpha = SX.sym("alpha")
k     = SX.sym("k")
sigma = SX.sym("sigma")
params = [eps,mu,alpha,k,sigma]

#! Variables
a     = SX.sym("a")
gamma = SX.sym("gamma")

#! Equations
res0 = mu*a+1.0/2*k*a*sin(gamma)
res1 = -sigma * a + 3.0/4*alpha*a**3+k*a*cos(gamma)

#! Numerical values
sigma_ = 0.1
alpha_ = 0.1
k_     = 0.2
params_ = [0.1,0.1,alpha_,k_,sigma_]

#! We create a NLPImplicitSolver instance
f=Function("f", [vertcat(a, gamma), vertcat(*params)], [vertcat(res0, res1)])
opts = {}
opts["nlpsol"] = "ipopt"
opts["nlpsol_options"] = {"ipopt.tol":1e-14}
s=rootfinder("s", "nlpsol", f, opts)

#$ Initialize [$a$,$\gamma$] with a guess and solve
x_ = s([1,-1], params_)
print("Solution = ", x_)

#! Compare with the analytic solution:
x = [sqrt(4.0/3*sigma_/alpha_),-0.5*pi]
print("Reference solution = ", x)

#! We show that the residual is indeed (close to) zero
residual = f(x_, params_)
print("residual = ", residual)

for i in range(1):
  assert(abs(x_[i]-x[i])<1e-6)

#! Solver statistics
print(s.stats())
