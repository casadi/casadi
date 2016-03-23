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
#! Integrator tolerances
#! =====================
from casadi import *
from numpy import *
from pylab import *

x=SX.sym('x') 
dx=SX.sym('dx')
states = vertcat(x,dx)

dae={'x':states, 'ode':vertcat(dx,-x)}

tend = 2*pi*3
ts = linspace(0,tend,1000)

tolerances = [-10,-5,-4,-3,-2,-1]

figure()

for tol in tolerances:
  opts = {'reltol':10.0**tol, 'abstol':10.0**tol, 'grid':ts, 'output_t0':True}
  F = integrator('F', 'cvodes', dae, opts)
  res = F(x0=[1,0])

  plot(ts,array(res['xf'])[0,:].T,label='tol = 1e%d' % tol)

legend( loc='upper left')
xlabel('Time [s]')
ylabel('State x [-]')
show()


tolerances = logspace(-15,1,500)
endresult=[]

for tol in tolerances:
  opts = {}
  opts['reltol'] = tol
  opts['abstol'] = tol
  opts['tf'] = tend
  F = integrator('F', 'cvodes', dae, opts)
  res = F(x0=[1,0])
  endresult.append(res['xf'][0])
  
figure()
loglog(tolerances,(array(endresult)-1),'b',label='Positive error')
loglog(tolerances,-(array(endresult)-1),'r',label='Negative error')
xlabel('Integrator relative tolerance')
ylabel('Error at the end of integration time')
legend(loc='upper left')
show()
