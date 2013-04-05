# -*- coding: utf-8 -*-
"""
Created on Fri May 11 16:46:10 2012

@author: janderss
"""

from casadi import *
x=ssym("x",3,1); u=ssym("u",1,1)
f =      [(1-x[1]*x[1])*x[0]-x[1]+u,\
          x[0],\
          x[0]*x[0]+x[1]*x[1]+u*u]
T      = 10.0;       E      = x[2]
x0     = [0,1,0];
xf_min = [0,0,-inf]; xf_max = [0,0,inf]
u_min  = -0.75;      u_max  = 1.0
ffcn = SXFunction([x,u],[f]); ffcn.init()
Efcn = SXFunction([x,u],[E]); Efcn.init()
N = 20; U = msym("U",N,1)
xk=msym("xk",3,1); uk=msym("uk",1,1)
M = 100; dt = T/(M*N); xkj = xk
for j in range(M):
  [f_kj] = ffcn.call([xkj,uk])
  xkj = xkj + dt*f_kj
integrator = MXFunction([xk,uk],[xkj])
integrator.init()
Xk = x0
for k in range(N):
  [Xk] = integrator.call([Xk,U[k]])
[obj] = Efcn.call([Xk,U[N-1]])
f_NLP = MXFunction([U],[obj])
con = vertcat([Xk[0],Xk[1]])
g_NLP = MXFunction([U],[con])
solver = IpoptSolver(f_NLP,g_NLP)
solver.setOption("max_iter",100)
#solver.setOption("expand_f",True)
#solver.setOption("expand_g",True)
solver.setOption("generate_hessian",True)
solver.setOption("derivative_test","second-order")
solver.init()
solver.setInput(N*[1],"x_init")
solver.setInput(N*[u_min],"lbx")
solver.setInput(N*[u_max],"ubx")
solver.setInput([0,0],"lbg")
solver.setInput([0,0],"ubg")
solver.solve()
