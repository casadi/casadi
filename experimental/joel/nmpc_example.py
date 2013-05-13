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
from casadi import *

# Formulate the NLP
N = 30
u = msym("u",N); p = msym("p")
f = inner_prod(u,u)
g = MX([]); gmax = []; gmin = []
x = p
for k in range(N):
  x = x+0.1*(x*(x+1)+u[k])
  x.lift(0.0)
  f = f + x*x
  g.append(x)
  gmin.append(-1 if k<N-1 else 0)
  gmax.append( 1 if k<N-1 else 0)

# Allocate NLP solver
h = MXFunction(nlpIn(x=u,p=p),nlpOut(f=f,g=g));
S = SCPgen(h)
S.setOption("qp_solver",QPOasesSolver)
S.setOption("qp_solver_options",{"printLevel":"none"}) # Should be automatic
S.init()

# Pass bounds and solve
S.setInput( 0.3,"p")    # S.setInput( 0.3,"p")
S.setInput(-1.0,"lbx")  # S.setInput(-1.0,"lbx")
S.setInput( 1.0,"ubx")  # S.setInput( 1.0,"ubx")
S.setInput(gmin,"lbg")  # S.setInput(gmin,"lbg")
S.setInput(gmax,"ubg")  # S.setInput(gmax,"ubg")
S.solve()

# Visualize the trajectory
from matplotlib.pylab import *
u = S.output("x")
plot(u)
x = DMatrix(0.30)
for k in range(30):
    xk = x[-1,0]
    x.append(xk + 0.1*(xk*(xk+1) + u[k]))
plot(x)
legend(['u','x'])
grid()
show()

