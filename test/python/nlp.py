from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *

class NLPtests(casadiTestCase):
  def testIPOPT(self):
    self.message("trivial IPOP")
    x=SX("x")
    f=SXFunction([x],[(x-1)**2])
    g=SXFunction([x],[x])
    
    solver = IpoptSolver(f,g)
    solver.setOption("tol",1e-5)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",100)
    solver.setOption("print_level",0)
    solver.init()
    solver.input(NLP_LBX).set([-10])
    solver.input(NLP_UBX).set([10])
    solver.input(NLP_LBG).set([-10])
    solver.input(NLP_UBG).set([10])
    solver.solve()
    self.assertAlmostEqual(solver.output(NLP_COST)[0],0,10,"IPOPT")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],1,10,"IPOPT")

  def testIPOPTmx(self):
    self.message("trivial IPOP, using MX")
    x=MX("x")
    f=MXFunction([x],[(x-1)**2])
    g=MXFunction([x],[2*x])
    
    solver = IpoptSolver(f,g)
    solver.setOption("tol",1e-10)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",100)
    solver.setOption("print_level",0)
    solver.init()
    solver.input(NLP_LBX).set([-10])
    solver.input(NLP_UBX).set([10])
    solver.input(NLP_LBG).set([-10])
    solver.input(NLP_UBG).set([10])
    solver.solve()
    self.assertAlmostEqual(solver.output(NLP_COST)[0],0,10,"IPOPT")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],1,10,"IPOPT")
    
  def testIPOPTc(self):
    self.message("trivial IPOP, overconstrained")
    x=SX("x")
    f=SXFunction([x],[(x-1)**2])
    g=SXFunction([x],[[x,x,x]])
    
    solver = IpoptSolver(f,g)
    solver.setOption("tol",1e-5)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",100)
    solver.setOption("print_level",0)
    solver.init()
    solver.input(NLP_LBX).set([-10])
    solver.input(NLP_UBX).set([10])
    solver.input(NLP_LBG).set([-10, -10, -10])
    solver.input(NLP_UBG).set([10, 10, 10])
    solver.solve()
    self.assertAlmostEqual(solver.output(NLP_COST)[0],0,10,"IPOPT")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],1,10,"IPOPT")
    
  def testIPOPTc2(self):
    self.message("trivial IPOP, overconstrained")
    x=SX("x")
    f=SXFunction([x],[(x-1)**2])
    g=SXFunction([x],[[x,x,x+x]])
    
    solver = IpoptSolver(f,g)
    solver.setOption("tol",1e-10)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",100)
    solver.setOption("print_level",0)
    solver.init()
    solver.input(NLP_LBX).set([-10])
    solver.input(NLP_UBX).set([10])
    solver.input(NLP_LBG).set([-10, -10, -10])
    solver.input(NLP_UBG).set([10, 10, 10])
    solver.solve()
    self.assertAlmostEqual(solver.output(NLP_COST)[0],0,10,"IPOPT")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],1,8,"IPOPT")
    
  def testIPOPTcmx(self):
    self.message("trivial IPOP, overconstrained, using MX")
    x=MX("x")
    f=MXFunction([x],[(x-1)**2])
    g=MXFunction([x],[vertcat([2*x,3*x,4*x])])
    
    solver = IpoptSolver(f,g)
    solver.setOption("tol",1e-10)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",100)
    solver.setOption("print_level",0)
    solver.init()
    solver.input(NLP_LBX).set([-10])
    solver.input(NLP_UBX).set([10])
    solver.input(NLP_LBG).set([-10,-10,-10])
    solver.input(NLP_UBG).set([10,10,10])
    solver.solve()
    self.assertAlmostEqual(solver.output(NLP_COST)[0],0,10,"IPOPT")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],1,10,"IPOPT")
    
if __name__ == '__main__':
    unittest.main()

