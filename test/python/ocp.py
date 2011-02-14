from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *

class OCPtests(casadiTestCase):
  def testdiscrete(self):
    self.message("Linear-quadratic problem, discrete, using IPOPT")
    # inspired by www.cs.umsl.edu/~janikow/publications/1992/GAforOpt/text.pdf
    a=1.0
    b=1.0
    q=1.0
    s=1.0
    r=1.0
    x0=100
    
    N=100
    
    X=symbolic("X",N+1)
    U=symbolic("U",N)
    
    V = vstack([X,U])
    
    cost = 0
    for i in range(N):
      cost = cost + s*X[i]**2+r*U[i]**2
    cost = cost + q*X[N]**2
    
    f = SXFunction([V],[cost])
    
    g = SXFunction([V],[vertcat([X[0]-x0,X[1:,0]-(a*X[:-1,0]+b*U)])])
    
    solver = IpoptSolver(f,g)
    solver.setOption("tol",1e-5)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",100)
    solver.setOption("print_level",0)
    solver.init()
    solver.input(NLP_LBX).set([-1000 for i in range(V.size)])
    solver.input(NLP_UBX).set([1000 for i in range(V.size)])
    solver.input(NLP_LBG).set([0 for i in range(N+1)])
    solver.input(NLP_UBG).set([0 for i in range(N+1)])
    solver.solve()
    ocp_sol=solver.output(NLP_COST)[0]
    # solve the ricatti equation exactly
    K = q+0.0
    for i in range(N):
      K = s+r*a**2*K/(r+b**2*K)
    exact_sol=K * x0**2
    self.assertAlmostEqual(ocp_sol,exact_sol,10,"Linear-quadratic problem solution using IPOPT")
    
    

if __name__ == '__main__':
    unittest.main()

