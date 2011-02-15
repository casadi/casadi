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
    f=SXFunction([x],[[(x-1)**2]])
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
    
    
  def testIPOPTrb(self):
    self.message("IPOPT rosenbrock, limited-memorey hessian approx")
    x=SX("x")
    y=SX("y")
    
    f=SXFunction([[x,y]],[(1-x)**2+100*(y-x**2)**2])
    
    solver = IpoptSolver(f)
    solver.setOption("tol",1e-8)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",100)
    #solver.setOption("print_level",0)
    solver.init()
    solver.input(NLP_LBX).set([-10]*2)
    solver.input(NLP_UBX).set([10]*2)
    solver.solve()
    self.assertAlmostEqual(solver.output(NLP_COST)[0],0,10,"IPOPT")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],1,10,"IPOPT")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[1],1,10,"IPOPT")

    
  def testIPOPTrb(self):
    self.message("IPOPT rosenbrock, limited-memorey hessian approx")
    x=SX("x")
    y=SX("y")
    
    f=SXFunction([[x,y]],[(1-x)**2+100*(y-x**2)**2])
    g=SXFunction([[x,y]],[x+y])
    solver = IpoptSolver(f,g)
    solver.setOption("tol",1e-8)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",100)
    #solver.setOption("print_level",0)
    solver.init()
    solver.input(NLP_LBX).set([-10]*2)
    solver.input(NLP_UBX).set([10]*2)
    solver.input(NLP_LBG).set([-10])
    solver.input(NLP_UBG).set([10])
    solver.solve()
    self.assertAlmostEqual(solver.output(NLP_COST)[0],0,10,"IPOPT")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],1,10,"IPOPT")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[1],1,10,"IPOPT")
    
    
  def testIPOPTrhb2(self):
    return
    self.message("IPOPT rosenbrock, exact hessian")
    x=SX("x")
    y=SX("y")
    
    f=SXFunction([[x,y]],[(1-x)**2+100*(y-x**2)**2])
    g=SXFunction([[x,y]],[x+y])
    solver = IpoptSolver(f,g)
    
    h=f.hessian()
    h.init()
    h.input().set([1,1])
    h.evaluate()
    print array(h.output())
    
    
    solver = IpoptSolver(f,g,h)
    solver.setOption("tol",1e-8)
    solver.setOption("max_iter",100)
    solver.setOption("hessian_approximation", "exact")
    #solver.setOption("print_level",0)
    solver.init()
    solver.input(NLP_LBX).set([-10]*2)
    solver.input(NLP_UBX).set([10]*2)
    solver.input(NLP_LBG).set([-10])
    solver.input(NLP_UBG).set([10])
    solver.solve()
    self.assertAlmostEqual(solver.output(NLP_COST)[0],0,10,"IPOPT")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],1,10,"IPOPT")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[1],1,10,"IPOPT")
     
  def testIPOPTrhb(self):
    self.message("IPOPT rosenbrock, exact hessian")
    x=SX("x")
    y=SX("y")
    
    obj=(1-x)**2+100*(y-x**2)**2
    f=SXFunction([[x,y]],[obj])
    
    sigma=SX("sigma")
    
    h=SXFunction([[x,y],[],[sigma]],[sigma*hessian(obj,[x,y])])
    
    solver = IpoptSolver(f,FX(),h)
    solver.setOption("tol",1e-10)
    solver.setOption("max_iter",100)
    solver.setOption("hessian_approximation", "exact")
    #solver.setOption("print_level",0)
    solver.init()
    solver.input(NLP_LBX).set([-10]*2)
    solver.input(NLP_UBX).set([10]*2)
    solver.solve()
    self.assertAlmostEqual(solver.output(NLP_COST)[0],0,10,"IPOPT")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],1,10,"IPOPT")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[1],1,10,"IPOPT")
     
  def testIPOPTnorm(self):
    return 
    self.message("IPOPT min ||x||_2")
    def norm_2(mx):
      return sqrt(inner_prod(mx,mx))
    N=10
    x=MX("x",N)
    x0=linspace(0,1,N)
    X0=MX(x0)
    f=MXFunction([x],[norm_2(x-X0)])
    g=MXFunction([x],[2*x])
    solver = IpoptSolver(f,g)
    solver.setOption("tol",1e-5)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",103)
    #solver.setOption("print_level",0)
    solver.init()
    solver.input(NLP_LBX).set([-10]*N)
    solver.input(NLP_UBX).set([10]*N)
    solver.input(NLP_LBG).set([-10]*N)
    solver.input(NLP_UBG).set([10]*N)
    solver.solve()
    print "residuals"
    print array(solver.output(NLP_X_OPT)).squeeze()-x0
    #self.assertAlmostEqual(solver.output(NLP_COST)[0],0,10,"IPOPT")
    #self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],1,10,"IPOPT")
    
  def testIPOPTnormSX(self):
    self.message("IPOPT min ||x||_2")
    def norm_2(mx):
      return inner_prod(mx,mx)
    N=10
    x=symbolic("x",N)
    x0=linspace(0,1,N)
    X0=DMatrix(x0)
    print norm_2(x-X0)
    f=SXFunction([x],[norm_2(x-X0)])
    #g=SXFunction([x],[2*x])
    solver = IpoptSolver(f)
    solver.setOption("tol",1e-5)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",103)
    solver.setOption("derivative_test","first-order")
    #solver.setOption("print_level",0)
    solver.init()
    solver.input(NLP_LBX).set([-10]*N)
    solver.input(NLP_UBX).set([10]*N)
    #solver.input(NLP_LBG).set([-10]*N)
    #solver.input(NLP_UBG).set([10]*N)
    solver.solve()
    print "residuals"
    print array(solver.output(NLP_X_OPT)).squeeze()-x0
    #self.assertAlmostEqual(solver.output(NLP_COST)[0],0,10,"IPOPT")
    #self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],1,10,"IPOPT")
  
  def testIPOPTnoc(self):
    self.message("trivial IPOP, no constraints")
    """ There is an assertion error thrown, but still it works"""
    x=SX("x")
    f=SXFunction([x],[(x-1)**2])
    solver = IpoptSolver(f)
    solver.setOption("tol",1e-10)
    solver.setOption("print_level",0)
    solver.init()
    solver.input(NLP_LBX).set([-10])
    solver.input(NLP_UBX).set([10])
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

  def testIPOPTdeg(self):
    self.message("degenerate optimization IPOP")
    x=SX("x")
    y=SX("y")
    f=SXFunction([[x,y]],[0])
    g=SXFunction([[x,y]],[[x-y,x]])
    
    solver = IpoptSolver(f,g)
    solver.setOption("tol",1e-5)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",100)
    solver.setOption("print_level",0)
    solver.init()
    solver.input(NLP_LBX).set([-10, -10])
    solver.input(NLP_UBX).set([10, 10])
    solver.input(NLP_LBG).set([0, 3])
    solver.input(NLP_UBG).set([0, 3])
    solver.solve()
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],solver.output(NLP_X_OPT)[1],10,"IPOPT")

  def testIPOPTdegc(self):
    self.message("degenerate optimization IPOP, overconstrained")
    x=SX("x")
    y=SX("y")
    f=SXFunction([[x,y]],[0])
    g=SXFunction([[x,y]],[[x-y,x,x+y]])
    
    solver = IpoptSolver(f,g)
    solver.setOption("tol",1e-5)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",100)
    solver.setOption("print_level",0)
    solver.init()
    solver.input(NLP_LBX).set([-10, -10])
    solver.input(NLP_UBX).set([10, 10])
    solver.input(NLP_LBG).set([0, 3 , -10])
    solver.input(NLP_UBG).set([0, 3, 10])
    solver.solve()
    # todo: catch error when set([0, 3 , 5]) two times
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],solver.output(NLP_X_OPT)[1],10,"IPOPT")
    
    
    
if __name__ == '__main__':
    unittest.main()

