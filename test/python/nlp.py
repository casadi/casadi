from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *

solvers= []
try:
  solvers.append(IpoptSolver)
except:
  pass
  
try:
  solvers.append(WorhpSolver)
except:
  pass

class NLPtests(casadiTestCase):
  def testIPOPT(self):
    
    x=SX("x")
    f=SXFunction([x],[[(x-1)**2]])
    g=SXFunction([x],[x])
    
    for Solver in solvers:
      self.message("trivial " + str(Solver))
      solver = Solver(f,g)
      for k,v in ({"tol":1e-5,"hessian_approximation":"limited-memory","max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
       
      solver.init()
      solver.input(NLP_LBX).set([-10])
      solver.input(NLP_UBX).set([10])
      solver.input(NLP_LBG).set([-10])
      solver.input(NLP_UBG).set([10])
      solver.solve()
      self.assertAlmostEqual(solver.output(NLP_COST)[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],1,10,str(Solver))
    
  def testIPOPTinf(self):
    self.message("trivial IPOPT, infinity bounds")
    x=SX("x")
    f=SXFunction([x],[[(x-1)**2]])
    g=SXFunction([x],[x])
    
    solver = IpoptSolver(f,g)
    solver.setOption("tol",1e-5)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",100)
    solver.setOption("print_level",0)
    solver.setOption("derivative_test","first-order")
    solver.init()
    solver.input(NLP_LBX).set([-Inf])
    solver.input(NLP_UBX).set([Inf])
    solver.input(NLP_LBG).set([-Inf])
    solver.input(NLP_UBG).set([Inf])
    solver.solve()
    self.assertAlmostEqual(solver.output(NLP_COST)[0],0,10,"IPOPT")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],1,10,"IPOPT")
    
  def testIPOPTrb(self):
    self.message("rosenbrock, limited-memorey hessian approx")
    x=SX("x")
    y=SX("y")
    
    f=SXFunction([[x,y]],[(1-x)**2+100*(y-x**2)**2])
    
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(f)
      for k,v in ({"tol":1e-8,"hessian_approximation":"limited-memory","max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.input(NLP_LBX).set([-10]*2)
      solver.input(NLP_UBX).set([10]*2)
      solver.solve()
      self.assertAlmostEqual(solver.output(NLP_COST)[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],1,10,str(Solver))
      self.assertAlmostEqual(solver.output(NLP_X_OPT)[1],1,10,str(Solver))

    
  def testIPOPTrb(self):
    self.message("rosenbrock, limited-memorey hessian approx")
    x=SX("x")
    y=SX("y")
    
    f=SXFunction([[x,y]],[(1-x)**2+100*(y-x**2)**2])
    g=SXFunction([[x,y]],[x+y])
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(f,g)
      for k,v in ({"tol":1e-8,"TolOpti":1e-20,"hessian_approximation":"limited-memory","max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.input(NLP_LBX).set([-10]*2)
      solver.input(NLP_UBX).set([10]*2)
      solver.input(NLP_LBG).set([-10])
      solver.input(NLP_UBG).set([10])
      solver.solve()
      self.assertAlmostEqual(solver.output(NLP_COST)[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],1,9,str(Solver))
      self.assertAlmostEqual(solver.output(NLP_X_OPT)[1],1,9,str(Solver))
    
    
  def testIPOPTrhb2(self):
    self.message("IPOPT rosenbrock, exact hessian, constrained")
    x=SX("x")
    y=SX("y")
    
    obj = (1-x)**2+100*(y-x**2)**2
    f=SXFunction([[x,y]],[obj])
    g=SXFunction([[x,y]],[x**2+y**2])
    solver = IpoptSolver(f,g)
    
    c_r = 4.56748075136258e-02;
    x_r = [7.86415156987791e-01,6.17698316967954e-01]
    
    sigma=SX("sigma")
    lambd=SX("lambd")
    h=SXFunction([[x,y],[lambd],[sigma]],[sigma*hessian(obj,[x,y])+lambd*hessian(g.outputSX(0),[x,y])])
    
    solver = IpoptSolver(f,g,h)
    solver.setOption("tol",1e-10)
    solver.setOption("max_iter",100)
    solver.setOption("hessian_approximation", "exact")
    solver.setOption("derivative_test","first-order")
    #solver.setOption("print_level",0)
    solver.init()
    solver.input(NLP_LBX).set([-10]*2)
    solver.input(NLP_UBX).set([10]*2)
    solver.input(NLP_LBG).set([0])
    solver.input(NLP_UBG).set([1])
    solver.solve()
    self.assertAlmostEqual(solver.output(NLP_COST)[0],c_r,10,"IPOPT")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],x_r[0],10,"IPOPT")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[1],x_r[1],10,"IPOPT")
    self.message(":warmstart")
    oldsolver=solver
    solver = IpoptSolver(f,g,h)
    solver.setOption("tol",1e-10)
    solver.setOption("max_iter",100)
    solver.setOption("hessian_approximation", "exact")
    #solver.setOption("print_level",0)
    solver.setOption("warm_start_init_point","yes")
    solver.setOption("warm_start_bound_push",1e-6)
    solver.setOption("warm_start_slack_bound_push",1e-6)
    solver.setOption("warm_start_mult_bound_push",1e-6)
    solver.setOption("mu_init",1e-6)
    solver.init()
    solver.input(NLP_LBX).set([-10]*2)
    solver.input(NLP_UBX).set([10]*2)
    solver.input(NLP_LBG).set([0])
    solver.input(NLP_UBG).set([1])
    solver.input(NLP_X_INIT).set(oldsolver.output(NLP_X_OPT))
    solver.input(NLP_LAMBDA_INIT).set(oldsolver.output(NLP_LAMBDA_G))
    solver.output(NLP_LAMBDA_X).set(oldsolver.output(NLP_LAMBDA_X))
    
    solver.solve()
    
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
    solver.setOption("derivative_test","first-order")
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
    self.message("IPOPT min ||x||^2_2")
    def norm_2(mx):
      return inner_prod(mx,mx)
    N=10
    x=MX("x",N)
    x0=linspace(0,1,N)
    X0=MX(x0)
    f=MXFunction([x],[norm_2(x-X0)])
    g=MXFunction([x],[2*x])
    solver = IpoptSolver(f,g)
    solver.setOption("tol",1e-5)
    #solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",103)
    solver.setOption("print_level",0)
    solver.setOption("derivative_test","first-order")
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
    self.message("IPOPT min ||x||^2_2")
    def norm_2(mx):
      return inner_prod(mx,mx)
    N=10
    x=ssym("x",N)
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
    solver.setOption("print_level",0)
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
    self.message("trivial IPOPT, no constraints")
    """ There is an assertion error thrown, but still it works"""
    x=ssym("x")
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
    self.message("trivial IPOPT, using MX")
    x=MX("x")
    f=MXFunction([x],[(x-1)**2])
    g=MXFunction([x],[2*x])
    
    solver = IpoptSolver(f,g)
    solver.setOption("tol",1e-10)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",100)
    solver.setOption("print_level",0)
    solver.setOption("derivative_test","first-order")
    solver.init()
    solver.input(NLP_LBX).set([-10])
    solver.input(NLP_UBX).set([10])
    solver.input(NLP_LBG).set([-10])
    solver.input(NLP_UBG).set([10])
    solver.solve()
    self.assertAlmostEqual(solver.output(NLP_COST)[0],0,10,"IPOPT")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],1,10,"IPOPT")
    
  def testIPOPTc(self):
    self.message("trivial IPOPT, overconstrained")
    x=SX("x")
    f=SXFunction([x],[(x-1)**2])
    g=SXFunction([x],[[x,x,x]])
    
    solver = IpoptSolver(f,g)
    solver.setOption("tol",1e-5)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",100)
    solver.setOption("print_level",0)
    solver.setOption("derivative_test","first-order")
    solver.init()
    solver.input(NLP_LBX).set([-10])
    solver.input(NLP_UBX).set([10])
    solver.input(NLP_LBG).set([-10, -10, -10])
    solver.input(NLP_UBG).set([10, 10, 10])
    solver.solve()
    self.assertAlmostEqual(solver.output(NLP_COST)[0],0,10,"IPOPT")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],1,10,"IPOPT")
    
  def testIPOPTc2(self):
    self.message("trivial IPOPT, overconstrained")
    x=SX("x")
    f=SXFunction([x],[(x-1)**2])
    g=SXFunction([x],[[x,x,x+x]])
    
    solver = IpoptSolver(f,g)
    solver.setOption("tol",1e-10)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",100)
    solver.setOption("print_level",0)
    solver.setOption("derivative_test","first-order")
    solver.init()
    solver.input(NLP_LBX).set([-10])
    solver.input(NLP_UBX).set([10])
    solver.input(NLP_LBG).set([-10, -10, -10])
    solver.input(NLP_UBG).set([10, 10, 10])
    solver.solve()
    self.assertAlmostEqual(solver.output(NLP_COST)[0],0,10,"IPOPT")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],1,8,"IPOPT")
    
  def testIPOPTcmx(self):
    self.message("trivial IPOPT, overconstrained, using MX")
    x=MX("x")
    f=MXFunction([x],[(x-1)**2])
    g=MXFunction([x],[vertcat([2*x,3*x,4*x])])
    
    solver = IpoptSolver(f,g)
    solver.setOption("tol",1e-10)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",100)
    solver.setOption("print_level",0)
    solver.setOption("derivative_test","first-order")
    solver.init()
    solver.input(NLP_LBX).set([-10])
    solver.input(NLP_UBX).set([10])
    solver.input(NLP_LBG).set([-10,-10,-10])
    solver.input(NLP_UBG).set([10,10,10])
    solver.solve()
    self.assertAlmostEqual(solver.output(NLP_COST)[0],0,10,"IPOPT")
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],1,10,"IPOPT")

  def testIPOPTdeg(self):
    self.message("degenerate optimization IPOPT")
    x=SX("x")
    y=SX("y")
    f=SXFunction([[x,y]],[0])
    g=SXFunction([[x,y]],[[x-y,x]])
    
    solver = IpoptSolver(f,g)
    solver.setOption("tol",1e-5)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",100)
    solver.setOption("print_level",0)
    solver.setOption("derivative_test","first-order")
    solver.init()
    solver.input(NLP_LBX).set([-10, -10])
    solver.input(NLP_UBX).set([10, 10])
    solver.input(NLP_LBG).set([0, 3])
    solver.input(NLP_UBG).set([0, 3])
    solver.solve()
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],solver.output(NLP_X_OPT)[1],10,"IPOPT")

  def testIPOPTdegc(self):
    self.message("degenerate optimization IPOPT, overconstrained")
    x=SX("x")
    y=SX("y")
    f=SXFunction([[x,y]],[0])
    g=SXFunction([[x,y]],[[x-y,x,x+y]])
    
    solver = IpoptSolver(f,g)
    solver.setOption("tol",1e-5)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",100)
    solver.setOption("print_level",0)
    solver.setOption("derivative_test","first-order")
    solver.init()
    solver.input(NLP_LBX).set([-10, -10])
    solver.input(NLP_UBX).set([10, 10])
    solver.input(NLP_LBG).set([0, 3 , -10])
    solver.input(NLP_UBG).set([0, 3, 10])
    solver.solve()
    # todo: catch error when set([0, 3 , 5]) two times
    self.assertAlmostEqual(solver.output(NLP_X_OPT)[0],solver.output(NLP_X_OPT)[1],10,"IPOPT")
    
  def testKINSol1(self):
    self.message("Scalar KINSol problem, n=0")
    x=SX("x")
    f=SXFunction([x],[sin(x)])
    f.init()
    solver=KinsolSolver(f,1)
    solver.init()
    solver.output().set(6)
    solver.solve()
    self.assertAlmostEqual(solver.output()[0],2*pi,5)

  def testKINSol2(self):
    self.message("Scalar KINSol problem, n=1")
    x=SX("x")
    y=SX("y")
    n=0.2
    f=SXFunction([y,x],[sin(x)-y])
    f.init()
    solver=KinsolSolver(f,1)
    solver.setOption("linear_solver_creator",CSparse) # NOTE by Joel: Sensitivities of an implicit function requires a user-provided linear solver 
    solver.init()
    solver.fwdSeed().set(1)
    solver.adjSeed().set(1)
    solver.input().set(n)
    solver.evaluate(1,1)
    self.assertAlmostEqual(solver.output()[0],sin(n),6)
    self.assertAlmostEqual(solver.fwdSens()[0],cos(n),6)
    self.assertAlmostEqual(solver.adjSens()[0],cos(n),6)

  def testKINSol1c(self):
    self.message("Scalar KINSol problem, n=0, constraint")
    x=SX("x")
    f=SXFunction([x],[sin(x)])
    f.init()
    solver=KinsolSolver(f,1)
    solver.setOption("constraints",[-1])
    print solver.dictionary()
    solver.init()
    solver.output().set(-6)
    solver.solve()
    self.assertAlmostEqual(solver.output()[0],-2*pi,5)
    
if __name__ == '__main__':
    unittest.main()

