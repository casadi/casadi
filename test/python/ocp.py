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
    
  def test_singleshooting(self):
    self.message("Single shooting")
    p0 = 0.2
    y0= 1
    yc0=dy0=0
    te=0.4

    t=symbolic("t")
    q=symbolic("y",2,1)
    p=symbolic("p",1,1)
    # y
    # y'
    f=SXFunction([t,q,p],[vertcat([q[1],p[0]+q[1]**2 ])])
    f.init()
    
    integrator = CVodesIntegrator(f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    integrator.setOption("verbose",True)
    integrator.setOption("steps_per_checkpoint",10000)

    integrator.init()

    var = MX("var",2,1)
    par = MX("par",1,1)
    parMX= par
    
    t0   = MX(0)
    tend = MX(te)
    q0   = vertcat([var[0],par])
    par  = var[1]
    qend=integrator([t0,tend,q0,par,MX(2,1),MX()])
    
    parc = MX(0)
    
    f = MXFunction([var,parMX],[qend[0]])
    f.init()
    fc = MXFunction([var],[-f([var,parc])])
    fc.init()
    solver = IpoptSolver(fc)
    solver.setOption("tol",1e-12)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",10)
    solver.setOption("derivative_test","first-order")
    solver.setOption("print_level",0)
    solver.init()
    solver.input(NLP_LBX).set([-1, -1])
    solver.input(NLP_UBX).set([1, 0.2])
    solver.solve()
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],1,8,"X_opt")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[1],0.2,8,"X_opt")
    
    self.assertAlmostEqual(solver.output(NLP_LAMBDA_UBX)[0],1,8,"Cost should be linear in y0")
    self.assertAlmostEqual(solver.output(NLP_LAMBDA_UBX)[1],(sqrt(p0)*(te*yc0**2-yc0+p0*te)*tan(arctan(yc0/sqrt(p0))+sqrt(p0)*te)+yc0**2)/(2*p0*yc0**2+2*p0**2),8,"Cost should be linear in y0")
    self.assertAlmostEqual(-solver.output(NLP_COST)[0],(2*y0-log(yc0**2/p0+1))/2-log(cos(arctan(yc0/sqrt(p0))+sqrt(p0)*te)),7,"Cost")
    self.assertAlmostEqual(solver.output(NLP_LAMBDA_LBX)[0],0,8,"Constraint is supposed to be unactive")
    self.assertAlmostEqual(solver.output(NLP_LAMBDA_LBX)[1],0,8,"Constraint is supposed to be unactive")
  
  def test_singleshooting2(self):
    self.message("Single shooting 2")
    p0 = 0.2
    y0= 0.2
    yc0=dy0=0.1
    te=0.4

    t=symbolic("t")
    q=symbolic("y",2,1)
    p=symbolic("p",1,1)
    # y
    # y'
    f=SXFunction([t,q,p],[vertcat([q[1],p[0]+q[1]**2 ])])
    f.init()
    
    integrator = CVodesIntegrator(f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    integrator.setOption("verbose",True)
    integrator.setOption("steps_per_checkpoint",10000)

    integrator.init()

    var = MX("var",2,1)
    par = MX("par",1,1)
    
    t0   = MX(0)
    tend = MX(te)
    q0   = vertcat([var[0],par])
    parl  = var[1]
    qend=integrator([t0,tend,q0,parl,MX(2,1),MX()])
    
    parc = MX(dy0)
    
    f = MXFunction([var,par],[qend[0]])
    f.init()
    fc = MXFunction([var],[-f([var,parc])])
    fc.init()
    
    g = MXFunction([var],[var[0]-var[1]])
    g.init()
    
    solver = IpoptSolver(fc,g)
    solver.setOption("tol",1e-12)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",10)
    solver.setOption("derivative_test","first-order")
    #solver.setOption("print_level",0)
    solver.init()
    solver.input(NLP_LBX).set([-1, -1])
    solver.input(NLP_UBX).set([1, 0.2])
    solver.input(NLP_LBG).set([-1])
    solver.input(NLP_UBG).set([0])
    solver.solve()

    self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],0.2,6,"X_opt")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[1],0.2,6,"X_opt")
    
    self.assertAlmostEqual(solver.output(NLP_LAMBDA_UBX)[0],0,8,"Constraint is supposed to be unactive")
    dfdp0 = (sqrt(p0)*(te*yc0**2-yc0+p0*te)*tan(arctan(yc0/sqrt(p0))+sqrt(p0)*te)+yc0**2)/(2*p0*yc0**2+2*p0**2)
    self.assertAlmostEqual(solver.output(NLP_LAMBDA_UBX)[1],1+dfdp0,8)
    self.assertAlmostEqual(solver.output(NLP_LAMBDA_OPT)[0],1,8)
    self.assertAlmostEqual(-solver.output(NLP_COST)[0],(2*y0-log(yc0**2/p0+1))/2-log(cos(arctan(yc0/sqrt(p0))+sqrt(p0)*te)),7,"Cost")
    self.assertAlmostEqual(solver.output(NLP_LAMBDA_LBX)[0],0,8,"Constraint is supposed to be unactive")
    self.assertAlmostEqual(solver.output(NLP_LAMBDA_LBX)[1],0,8,"Constraint is supposed to be unactive")  
if __name__ == '__main__':
    unittest.main()

