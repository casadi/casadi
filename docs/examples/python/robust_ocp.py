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
from pylab import *

from casadi import *
from casadi.tools import *

try:
  DpleSolver.loadPlugin("slicot")
except:
  print "Slicot is not available as plugin"
  import sys
  sys.exit(0)

#  ==================================
#   Define the model xdot = f(x,u,w)
#  ==================================
xy = struct(["x","y"])
x  = struct_symSX([entry("p",struct=xy), entry("v",struct=xy)])
u  = struct_symSX(xy)
w  = struct_symSX(xy)

rhs = struct_SX(x)
rhs["p"]  = x["v"]
rhs["v"]  = -10*(x["p"]-u) - x["v"]*sqrt(sum_square(x["v"])+1) + w

T = SX.sym("T") # Time-period

f = SXFunction(daeIn(x=x,p=vertcat([T,u,w])),daeOut(ode=T*rhs))
f.init()

#  ==================================
#    Construct the path constraints
#  ==================================

gamma = 1

hyper = [ ( vertcat([1,1]),   vertcat([0,0  ]) , 4),
          ( vertcat([0.5,2]), vertcat([1,0.5]) , 4)]

P = SX.sym("P",4,4)

h_nom    = [ sumAll(((x["p"]-p)/s)**n) - 1 for s,p,n in hyper ]
h_margin = [ sqrt(quad_form(jacobian(i,x).T,P)) for i in h_nom ]

h_robust = SXFunction([x,P],[ n - gamma*m for n,m in zip(h_nom, h_margin) ])
h_robust.init()


#  ==================================
#  Construct an integrating block Phi
#  ==================================

N = 20
Tref = 4.0

ts = [i*1.0/N for i in range(N)]

Phi = Integrator("rk",f)
Phi.setOption("number_of_finite_elements",2)
Phi.setOption("tf",1.0/N)
Phi.init()

# Obtain some sensitivity functions from the integrating block
APhi = Phi.jacobian("x0","xf")
APhi.init()

CPhi = Phi.jacobian("p","xf")
CPhi.init()

#  ==================================
#         Construct the NLP
#  ==================================

V = struct_symMX([
      entry("X",struct=x,repeat=N+1),
      entry("U",struct=u,repeat=N),
      "T"
  ])

g_coupling = []
As = []
Qs = []
for k in range(N):
  pk = vertcat([V["T"],V["U",k],DMatrix.zeros(w.shape)])
  xkp = Phi(x0=V["X",k],p=pk)["xf"]
  
  g_coupling.append(xkp-V["X",k+1])
  
  A = APhi(x0=V["X",k],p=pk)["jac"]
  As.append(A)
  
  C = CPhi(x0=V["X",k],p=pk)["jac"][:,-2:]
  Qs.append(mul(C,C.T)*50)
 

# DPLE solver
dple = DpleSolver("slicot",[i.sparsity() for i in As],[i.sparsity() for i in Qs])
dple.setOption("linear_solver","csparse")
dple.init()

Ps = horzsplit(dple(a=horzcat(As),v=horzcat(Qs))["p"],x.shape[0])

# Constraints
g = struct_MX([
     entry("periodic", expr=V["X",0]-V["X",-1]),
     entry("coupling", expr=g_coupling),
     entry("obstacle", expr=[ h_robust([V["X",k],Ps[k]]) for k in range(N) ]),
  ])
 
# Objective function
F = V["T"] + 1e-2*sum_square(vertcat(V["U"]))/2/N

nlp = MXFunction(nlpIn(x=V),nlpOut(f=F,g=g))

#  ==================================
#   Set NLP initial guess and bounds
#  ==================================

solver = NlpSolver("ipopt",nlp)
solver.setOption("hessian_approximation","limited-memory")
solver.init()

V0 = V(0.0)

V0["T"] = Tref

for k in range(N+1):
  V0["X",k,"p","x"] = 3*sin(2*pi*k/N)
  V0["X",k,"p","y"] = 3*cos(2*pi*k/N)

solver.setInput(V0,"x0")

lbg = g(-inf);    ubg = g(inf)
lbg["coupling"] = ubg["coupling"] = 0
lbg["periodic"] = ubg["periodic"] = 0
lbg["obstacle"] = 0

solver.setInput(lbg,"lbg")
solver.setInput(ubg,"ubg")

#  ==================================
#              Solve NLP
#  ==================================

solver.evaluate()

sol = V(solver.getOutput())

#  ==================================
#              Plotting
#  ==================================

figure()

plot(sol["X",:,"p","x"],sol["X",:,"p","y"])

theta = linspace(0,2*pi,1000)
for s,p,n in hyper:
  fill(s[0]*fabs(cos(theta))**(2.0/n)*sign(cos(theta)) + p[0],s[1]*fabs(sin(theta))**(2.0/n)*sign(sin(theta)) + p[1],'r')

# Plot the unertainty ellipsoids
Pf = MXFunction(nlpIn(x=V),Ps)
Pf.init()

Pf.setInput(sol)
Pf.evaluate()

circle = array([[sin(x),cos(x)] for x in linspace(0,2*pi,100)]).T

for k in range(N):
  w,v = numpy.linalg.eig(Pf.getOutput(k)[:2,:2])
  W   = mul(v,diag(sqrt(w)))
  e = mul(W,circle)
  plot(e[0,:].T + sol["X",k,"p","x"],e[1,:].T + sol["X",k,"p","y"],'r')

show()

