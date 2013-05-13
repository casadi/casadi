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
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *

solvers= []
  
try:
  solvers.append(WorhpSolver)
  print "Will test WorhpSolver"
except:
  pass
  
try:
  solvers.append(IpoptSolver)
  print "Will test IpoptSolver"
except:
  pass

#try:
#  solvers.append(KnitroSolver)
#  print "Will test KnitroSolver"
#except:
#  pass
  
solvers.append(SQPMethod)
print "Will test SQPMethod"

qpsolver = NLPQPSolver
try:
  qpsolver_options = {"nlp_solver": IpoptSolver, "nlp_solver_options": {"tol": 1e-12} }
except:
  qpsolver_options = {}
#qpsolver = QPOasesSolver

class NLPtests(casadiTestCase):
  def testIPOPT(self):
    x=SX("x")
    nlp=SXFunction(nlIn(x=x),nlOut(f=(x-1)**2,g=x))
    
    for Solver in solvers:
      self.message("trivial " + str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-5,"hessian_approximation":"limited-memory","max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver, "qp_solver_options" : qpsolver_options }).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
       
      solver.init()
      solver.input("lbx").set([-10])
      solver.input("ubx").set([10])
      solver.input("lbg").set([-10])
      solver.input("ubg").set([10])
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,9,str(Solver))
      self.assertAlmostEqual(solver.output("g")[0],1,9,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[0],0,9,str(Solver))
      self.assertAlmostEqual(solver.output("lam_g")[0],0,9,str(Solver))
      
  def testIPOPT_par(self):
    x=SX("x")
    p=SX("p")
    nlp=SXFunction(nlIn(x=x,p=p),nlOut(f=(x-p)**2,g=x))
    
    for Solver in solvers:
      self.message("trivial " + str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-5,"hessian_approximation":"limited-memory","max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver, "qp_solver_options" : qpsolver_options}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.input("lbx").set([-10])
      solver.input("ubx").set([10])
      solver.input("lbg").set([-10])
      solver.input("ubg").set([10])
      solver.input("p").set(1)
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,9,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[0],0,9,str(Solver))
      self.assertAlmostEqual(solver.output("lam_g")[0],0,9,str(Solver))
      
  def testIPOPTinf(self):
    self.message("trivial IPOPT, infinity bounds")
    x=SX("x")
    nlp=SXFunction(nlIn(x=x),nlOut(f=(x-1)**2,g=x))
    
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-5,"hessian_approximation":"limited-memory","max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver, "qp_solver_options" : qpsolver_options}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.input("lbx").set([-Inf])
      solver.input("ubx").set([Inf])
      solver.input("lbg").set([-Inf])
      solver.input("ubg").set([Inf])
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,7,str(Solver) + str(solver.output("x")[0]-1))
      self.assertAlmostEqual(solver.output("lam_x")[0],0,9,str(Solver))
      self.assertAlmostEqual(solver.output("lam_g")[0],0,9,str(Solver))
      
  def testIPOPTrb(self):
    self.message("rosenbrock, limited-memory hessian approx")
    x=SX("x")
    y=SX("y")
    
    nlp=SXFunction(nlIn(x=vertcat([x,y])),nlOut(f=(1-x)**2+100*(y-x**2)**2))
    
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-9,"TolOpti":1e-14,"hessian_approximation":"limited-memory","maxiter":100,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver,"qp_solver_options" : qpsolver_options}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.input("lbx").set([-10]*2)
      solver.input("ubx").set([10]*2)
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,7,str(Solver))
      self.assertAlmostEqual(solver.output("x")[1],1,7,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[0],0,9,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[1],0,9,str(Solver))
    
  def testIPOPTrb2(self):
    self.message("rosenbrock, limited-memory hessian approx")
    x=SX("x")
    y=SX("y")
    
    nlp=SXFunction(nlIn(x=vertcat([x,y])),nlOut(f=(1-x)**2+100*(y-x**2)**2,g=x+y))
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-8,"TolOpti":1e-20,"hessian_approximation":"limited-memory","max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver,"qp_solver_options" : qpsolver_options, "maxiter": 1000}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.input("lbx").set([-10]*2)
      solver.input("ubx").set([10]*2)
      solver.input("lbg").set([-10])
      solver.input("ubg").set([10])
      solver.solve()
      
      digits = 6

      self.assertAlmostEqual(solver.output("f")[0],0,digits,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,digits,str(Solver))
      self.assertAlmostEqual(solver.output("x")[1],1,digits,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[0],0,9,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[1],0,9,str(Solver))
      self.assertAlmostEqual(solver.output("lam_g")[0],0,9,str(Solver))
      
  def testIPOPTrbf(self):
    self.message("rosenbrock fixed, limited-memory hessian approx")
    x=SX("x")
    y=SX("y")
    
    nlp=SXFunction(nlIn(x=vertcat([x,y])),nlOut(f=(1-x)**2+100*(y-x**2)**2,g=x+y))
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-8,"TolOpti":1e-20,"hessian_approximation":"limited-memory","max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver,"qp_solver_options" : qpsolver_options}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.input("x0").set([0,1])
      solver.input("lbx").set([-10,1])
      solver.input("ubx").set([10,1])
      solver.input("lbg").set([-10])
      solver.input("ubg").set([10])
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,7,str(Solver))
      self.assertAlmostEqual(solver.output("x")[1],1,7,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[0],0,6,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[1],0,6,str(Solver))
      self.assertAlmostEqual(solver.output("lam_g")[0],0,6,str(Solver))
      
  def testIPOPTrhb2(self):
    self.message("rosenbrock, exact hessian, constrained")
    x=SX("x")
    y=SX("y")
    
    obj = (1-x)**2+100*(y-x**2)**2
    nlp=SXFunction(nlIn(x=vertcat([x,y])),nlOut(f=obj,g=x**2+y**2))
    
    c_r = 4.56748075136258e-02;
    x_r = [7.86415156987791e-01,6.17698316967954e-01]
    
    sigma=SX("sigma")
    lambd=SX("lambd")
    h=SXFunction(hessLagIn(x=vertcat([x,y]),lam_f=sigma,lam_g=lambd),
                 hessLagOut(hess=sigma*hessian(obj,vertcat([x,y]))+lambd*hessian(nlp.outputExpr("g"),vertcat([x,y]))))
    h.init()
    h.input().set([0.5,0.5])
    h.input(1).set(-40)
    h.input(2).set(1)
    h.evaluate()
    print h.output()
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption("hess_lag",h)
      for k,v in ({"tol":1e-10,"TolOpti":1e-20,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "MaxIter": 100,"derivative_test":"second-order","qp_solver": qpsolver,"qp_solver_options" : qpsolver_options}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
          
      solver.init()
      solver.input("x0").set([0.5,0.5])
      solver.input("lbx").set([-10]*2)
      solver.input("ubx").set([10]*2)
      solver.input("lbg").set([0])
      solver.input("ubg").set([1])
      solver.solve()
      
      digits = 5
        
      self.assertAlmostEqual(solver.output("f")[0],c_r,digits,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],x_r[0],digits,str(Solver))
      self.assertAlmostEqual(solver.output("x")[1],x_r[1],digits,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[0],0,8,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[1],0,8,str(Solver))
      self.assertAlmostEqual(solver.output("lam_g")[0],0.12149655447670,6,str(Solver))
      
    self.message(":warmstart")
    oldsolver=solver
    solver = IpoptSolver(nlp)
    solver.setOption("hess_lag",h)
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
    solver.input("lbx").set([-10]*2)
    solver.input("ubx").set([10]*2)
    solver.input("lbg").set([0])
    solver.input("ubg").set([1])
    solver.input("x0").set(oldsolver.output("x"))
    solver.input("lam_g0").set(oldsolver.output("lam_g"))
    solver.output("lam_x").set(oldsolver.output("lam_x"))
    
    
    solver.solve()

  def testIPOPTrhb2_gen(self):
    self.message("rosenbrock, exact hessian generated, constrained")
    x=SX("x")
    y=SX("y")
    
    obj = (1-x)**2+100*(y-x**2)**2
    nlp=SXFunction(nlIn(x=vertcat([x,y])),nlOut(f=obj,g=x**2+y**2))
    
    c_r = 4.56748075136258e-02;
    x_r = [7.86415156987791e-01,6.17698316967954e-01]
    
    sigma=SX("sigma")
    lambd=SX("lambd")
  
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-12,"TolOpti":1e-20,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "MaxIter": 100,"print_level":1,"derivative_test":"second-order","qp_solver": qpsolver,"qp_solver_options" : qpsolver_options, "toldx": 1e-15, "tolgl": 1e-15, "maxiter" : 200}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
          
      solver.init()
      solver.input("x0").set([0.5,0.5])
      solver.input("lbx").set([-10]*2)
      solver.input("ubx").set([10]*2)
      solver.input("lbg").set([0])
      solver.input("ubg").set([1])
      solver.solve()
      
      digits = 5
      
      self.assertAlmostEqual(solver.output("f")[0],c_r,digits,str(Solver) + str(solver.output("f")[0]) + ":" + str(c_r))
      self.assertAlmostEqual(solver.output("x")[0],x_r[0],digits,str(Solver))
      self.assertAlmostEqual(solver.output("x")[1],x_r[1],digits,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[0],0,8,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[1],0,8,str(Solver))
      self.assertAlmostEqual(solver.output("lam_g")[0],0.12149655447670,6,str(Solver))
      
  def testIPOPTrhb2_par(self):
    self.message("rosenbrock, exact hessian, constrained, ")
    x=SX("x")
    y=SX("y")
    p=SX("p")
    
    obj = (p-x)**2+100*(y-x**2)**2
    nlp=SXFunction(nlIn(x=vertcat([x,y]),p=p),nlOut(f=obj,g=x**2+y**2))
    
    c_r = 4.56748075136258e-02;
    x_r = [7.86415156987791e-01,6.17698316967954e-01]
    
    sigma=SX("sigma")
    lambd=SX("lambd")
    h=SXFunction(hessLagIn(x=vertcat([x,y]),lam_f=sigma,lam_g=lambd),
                 hessLagOut(hess=sigma*hessian(obj,vertcat([x,y]))+lambd*hessian(nlp.outputExpr("g"),vertcat([x,y]))))

    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption("hess_lag",h)
      for k,v in ({"tol":1e-10,"TolOpti":1e-20,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "MaxIter": 100,"print_level":1,"derivative_test":"second-order","qp_solver": qpsolver,"qp_solver_options" : qpsolver_options}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.input("x0").set([0.5,0.5])
      solver.input("lbx").set([-10]*2)
      solver.input("ubx").set([10]*2)
      solver.input("lbg").set([0])
      solver.input("ubg").set([1])
      solver.input("p").set([1])
      solver.solve()
      
      digits = 5
        
      self.assertAlmostEqual(solver.output("f")[0],c_r,digits,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],x_r[0],digits,str(Solver))
      self.assertAlmostEqual(solver.output("x")[1],x_r[1],digits,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[0],0,8,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[1],0,8,str(Solver))
      self.assertAlmostEqual(solver.output("lam_g")[0],0.12149655447670,6,str(Solver))

  def testIPOPTrhb2_gen_par(self):
    self.message("rosenbrock, exact hessian generated, constrained, parametric")
    x=SX("x")
    y=SX("y")
    p=SX("p")
    
    obj = (p-x)**2+100*(y-x**2)**2
    nlp=SXFunction(nlIn(x=vertcat([x,y]),p=p),nlOut(f=obj,g=x**2+y**2))
    
    c_r = 4.56748075136258e-02;
    x_r = [7.86415156987791e-01,6.17698316967954e-01]
    
    sigma=SX("sigma")
    lambd=SX("lambd")
  
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-10,"TolOpti":1e-20,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "MaxIter": 100,"print_level":1,"derivative_test":"second-order","qp_solver": qpsolver,"qp_solver_options" : qpsolver_options}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
          
      solver.init()
      solver.input("x0").set([0.5,0.5])
      solver.input("lbx").set([-10]*2)
      solver.input("ubx").set([10]*2)
      solver.input("lbg").set([0])
      solver.input("ubg").set([1])
      solver.input("p").set([1])
      solver.solve()
      
      digits = 5

      self.assertAlmostEqual(solver.output("f")[0],c_r,digits,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],x_r[0],digits,str(Solver))
      self.assertAlmostEqual(solver.output("x")[1],x_r[1],digits,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[0],0,8,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[1],0,8,str(Solver))
      self.assertAlmostEqual(solver.output("lam_g")[0],0.12149655447670,6,str(Solver))
      
  def testIPOPTrhb(self):
    self.message("rosenbrock, exact hessian")
    x=SX("x")
    y=SX("y")
    
    obj=(1-x)**2+100*(y-x**2)**2
    nlp=SXFunction(nlIn(x=vertcat([x,y])),nlOut(f=obj))
    
    sigma=SX("sigma")
    
    h=SXFunction(hessLagIn(x=vertcat([x,y]),lam_f=sigma),
                 hessLagOut(hess=sigma*hessian(obj,vertcat([x,y]))))
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption("hess_lag",h)
      for k,v in ({"tol":1e-10,"TolOpti":1e-20,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver,"qp_solver_options" : qpsolver_options}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      #solver.setOption("verbose",True)
      solver.init()
      solver.input("lbx").set([-10]*2)
      solver.input("ubx").set([10]*2)
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,9,str(Solver))
      self.assertAlmostEqual(solver.output("x")[1],1,9,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[0],0,8,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[1],0,8,str(Solver))

  def testIPOPTrhb_gen(self):
    self.message("rosenbrock, exact hessian generated")
    x=SX("x")
    y=SX("y")
    
    obj=(1-x)**2+100*(y-x**2)**2
    nlp=SXFunction(nlIn(x=vertcat([x,y])),nlOut(f=obj))
    
    sigma=SX("sigma")
    
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-10,"TolOpti":1e-20,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver,"qp_solver_options" : qpsolver_options}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      #solver.setOption("verbose",True)
      solver.init()
      solver.input("lbx").set([-10]*2)
      solver.input("ubx").set([10]*2)
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,9,str(Solver))
      self.assertAlmostEqual(solver.output("x")[1],1,9,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[0],0,8,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[1],0,8,str(Solver))

  def testIPOPTrhb_gen_xnonfree(self):
    self.message("rosenbrock, exact hessian generated, non-free x")
    x=SX("x")
    y=SX("y")
    
    obj=(1-x)**2+100*(y-x**2)**2
    nlp=SXFunction(nlIn(x=vertcat([x,y])),nlOut(f=obj))
    
    sigma=SX("sigma")
    
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-10,"TolOpti":1e-20,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver,"qp_solver_options" : qpsolver_options}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      #solver.setOption("verbose",True)
      solver.init()
      solver.input("lbx").set([1,-10])
      solver.input("ubx").set([1,10])
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,9,str(Solver))
      self.assertAlmostEqual(solver.output("x")[1],1,6,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[0],0,6,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[1],0,6,str(Solver))
      
  def testIPOPTrhb_par(self):
    self.message("rosenbrock, exact hessian, parametric")
    x=SX("x")
    y=SX("y")
    
    p=SX("p")
    obj=(p-x)**2+100*(y-x**2)**2
    nlp=SXFunction(nlIn(x=vertcat([x,y]),p=p),nlOut(f=obj))
    
    sigma=SX("sigma")
    
    h=SXFunction(hessLagIn(x=vertcat([x,y]),p=p,lam_f=sigma),
                 hessLagOut(hess=sigma*hessian(obj,vertcat([x,y]))))
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption("hess_lag",h)
      for k,v in ({"tol":1e-10,"TolOpti":1e-20,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver,"qp_solver_options" : qpsolver_options}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      #solver.setOption("verbose",True)
      solver.init()
      solver.input("lbx").set([-10]*2)
      solver.input("ubx").set([10]*2)
      solver.input("p").set(1)
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,9,str(Solver))
      self.assertAlmostEqual(solver.output("x")[1],1,9,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[0],0,8,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[1],0,8,str(Solver))

  def testIPOPTrhb_gen_par(self):
    self.message("rosenbrock, exact hessian generated, parametric")
    x=SX("x")
    y=SX("y")
    
    p=SX("p")
    obj=(p-x)**2+100*(y-x**2)**2
    nlp=SXFunction(nlIn(x=vertcat([x,y]),p=p),nlOut(f=obj))
    
    sigma=SX("sigma")
    
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-10,"TolOpti":1e-20,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver,"qp_solver_options" : qpsolver_options}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      #solver.setOption("verbose",True)
      solver.init()
      solver.input("lbx").set([-10]*2)
      solver.input("ubx").set([10]*2)
      solver.input("p").set(1)
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,9,str(Solver))
      self.assertAlmostEqual(solver.output("x")[1],1,9,str(Solver))
      
  def testIPOPTnorm(self):
    self.message("IPOPT min ||x||^2_2")
    def norm_2(mx):
      return inner_prod(mx,mx)
    N=10
    x=msym("x",N)
    x0=linspace(0,1,N)
    X0=MX(x0)
    nlp=MXFunction(nlIn(x=x),nlOut(f=norm_2(x-X0),g=2*x))
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-8,"max_iter":103, "MaxIter": 103,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver,"qp_solver_options" : qpsolver_options}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.input("lbx").set([-10]*N)
      solver.input("ubx").set([10]*N)
      solver.input("lbg").set([-10]*N)
      solver.input("ubg").set([10]*N)
      solver.solve()
      print "residuals"
      print array(solver.output("x")).squeeze()-x0
      self.assertAlmostEqual(solver.output("f")[0],0,10,str(Solver))
      self.checkarray(array(solver.output("x")).squeeze(),x0,str(Solver))
      self.checkarray(solver.output("lam_x"),DMatrix([0]*10),8,str(Solver),digits=8)
      self.assertAlmostEqual(solver.output("lam_g")[1],0,8,str(Solver))
      
  def testIPOPTnoc(self):
    self.message("trivial IPOPT, no constraints")
    """ There is an assertion error thrown, but still it works"""
    x=ssym("x")
    nlp=SXFunction(nlIn(x=x),nlOut(f=(x-1)**2))
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-10,"max_iter":103, "MaxIter": 103,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver,"qp_solver_options" : qpsolver_options}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver = IpoptSolver(nlp)
      solver.init()
      solver.input("lbx").set([-10])
      solver.input("ubx").set([10])
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,9,str(Solver))
    
  def testIPOPTmx(self):
    self.message("trivial IPOPT, using MX")
    x=MX("x")
    nlp=MXFunction(nlIn(x=x),nlOut(f=(x-1)**2,g=2*x))
    
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-10,"max_iter":103, "MaxIter": 103,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver, "qp_solver_options" : qpsolver_options}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.input("lbx").set([-10])
      solver.input("ubx").set([10])
      solver.input("lbg").set([-10])
      solver.input("ubg").set([10])
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,9,str(Solver))
    
  def testIPOPTc(self):
    self.message("trivial, overconstrained")
    x=SX("x")
    nlp=SXFunction(nlIn(x=x),nlOut(f=(x-1)**2,g=vertcat([x,x,x])))
    
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-5,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver, "qp_solver_options" : qpsolver_options}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.input("lbx").set([-10])
      solver.input("ubx").set([10])
      solver.input("lbg").set([-10, -10, -10])
      solver.input("ubg").set([10, 10, 10])
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],0,9,str(Solver) )
      self.assertAlmostEqual(solver.output("x")[0],1,9,str(Solver))
    
  def testIPOPTc2(self):
    self.message("trivial2, overconstrained")
    x=SX("x")
    nlp=SXFunction(nlIn(x=x),nlOut(f=(x-1)**2,g=vertcat([x,x,x+x])))
    
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-10,"max_iter":100, "hessian_approximation": "limited-memory", "MaxIter": 100,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver, "qp_solver_options" : qpsolver_options}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.input("lbx").set([-10])
      solver.input("ubx").set([10])
      solver.input("lbg").set([-10, -10, -10])
      solver.input("ubg").set([10, 10, 10])
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,8,str(Solver))
    
  def testIPOPTcmx(self):
    self.message("trivial , overconstrained, using MX")
    x=MX("x")
    nlp=MXFunction(nlIn(x=x),nlOut(f=(x-1)**2,g=vertcat([2*x,3*x,4*x])))
    
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-10,"max_iter":100, "hessian_approximation": "limited-memory", "MaxIter": 100,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver, "qp_solver_options" : qpsolver_options}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.input("lbx").set([-10])
      solver.input("ubx").set([10])
      solver.input("lbg").set([-10,-10,-10])
      solver.input("ubg").set([10,10,10])
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],0,9,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,8,str(Solver))

  def testIPOPTdeg(self):
    self.message("degenerate optimization IPOPT")
    x=SX("x")
    y=SX("y")
    nlp=SXFunction(nlIn(x=vertcat([x,y])),nlOut(f=0,g=vertcat([x-y,x])))
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-5,"max_iter":100, "hessian_approximation": "limited-memory", "MaxIter": 100,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver,"qp_solver_options" : qpsolver_options}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.input("lbx").set([-10, -10])
      solver.input("ubx").set([10, 10])
      solver.input("lbg").set([0, 3])
      solver.input("ubg").set([0, 3])
      solver.solve()
      self.assertAlmostEqual(solver.output("x")[0],solver.output("x")[1],10,"IPOPT")

  def testIPOPTdegc(self):
    self.message("degenerate optimization IPOPT, overconstrained")
    x=SX("x")
    y=SX("y")
    nlp=SXFunction(nlIn(x=vertcat([x,y])),nlOut(f=0,g=vertcat([x-y,x,x+y])))
    
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-5,"max_iter":100, "hessian_approximation": "limited-memory", "MaxIter": 100,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver,"qp_solver_options" : qpsolver_options}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)

      solver.init()
      solver.input("lbx").set([-10, -10])
      solver.input("ubx").set([10, 10])
      solver.input("lbg").set([0, 3 , -10])
      solver.input("ubg").set([0, 3, 10])
      solver.solve()
      # todo: catch error when set([0, 3 , 5]) two times
      self.assertAlmostEqual(solver.output("x")[0],solver.output("x")[1],10,"IPOPT")
      
  def testXfreeChange(self):
    self.message("Change in X settings")
    x=SX("x")
    y=SX("y")
    
    nlp=SXFunction(nlIn(x=vertcat([x,y])),nlOut(f=(1-x)**2+100*(y-x**2)**2,g=x+y))
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-8,"TolOpti":1e-20,"hessian_approximation":"limited-memory","max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver,"qp_solver_options" : qpsolver_options}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.input("x0").set([0,1])
      solver.input("lbx").set([-10,-10])
      solver.input("ubx").set([10,10])
      solver.input("lbg").set([-10])
      solver.input("ubg").set([10])
      solver.solve()
      solver.input("lbx").set([-10,1])
      solver.input("ubx").set([10,1])
      solver.input("lbg").set([-10])
      solver.input("ubg").set([10])
      solver.solve()
      
      self.assertAlmostEqual(solver.output("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,7,str(Solver))
      self.assertAlmostEqual(solver.output("x")[1],1,7,str(Solver))

  def testactiveLBX(self):
    self.message("active LBX")
    x=SX("x")
    y=SX("y")
    
    nlp=SXFunction(nlIn(x=vertcat([x,y])),nlOut(f=(1-x)**2+100*(y-x**2)**2,g=x+y))
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-8,"TolOpti":1e-20,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver,"qp_solver_options" : qpsolver_options, "hessian_approximation": "exact", "UserHM": True}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.input("x0").set([0,1])
      solver.input("lbx").set([-10,1.2])
      solver.input("ubx").set([10,2])
      solver.input("lbg").set([-10])
      solver.input("ubg").set([10])
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],9.0908263002590e-3,6,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1.0952466252248,6,str(Solver))
      self.assertAlmostEqual(solver.output("x")[1],1.2,5,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[0],0,8,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[1],-8.6963632695079e-2,4,str(Solver))
      self.assertAlmostEqual(solver.output("lam_g")[0],0,8,str(Solver))

  def testactiveLBG(self):
    self.message("active LBG")
    x=SX("x")
    y=SX("y")
    
    nlp=SXFunction(nlIn(x=vertcat([x,y])),nlOut(f=(1-x)**2+100*(y-x**2)**2,g=x+y))
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-8,"TolOpti":1e-20,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver,"qp_solver_options" : qpsolver_options, "hessian_approximation": "exact", "UserHM": True }).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.input("x0").set([0,1])
      solver.input("lbx").set([-10,-10])
      solver.input("ubx").set([10,10])
      solver.input("lbg").set([2.2])
      solver.input("ubg").set([10])
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],4.252906468284e-3,6,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1.065181061847138,6,str(Solver))
      self.assertAlmostEqual(solver.output("x")[1],1.1348189166291160,6,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[0],0,8,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[1],0,4,str(Solver))
      self.assertAlmostEqual(solver.output("lam_g")[0],-4.1644422845712e-2,3,str(Solver))

  def testactiveUBG(self):
    self.message("active UBG")
    x=SX("x")
    y=SX("y")
    
    nlp=SXFunction(nlIn(x=vertcat([x,y])),nlOut(f=(1-x)**2+100*(y-x**2)**2,g=x+y))
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-8,"TolOpti":1e-20,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver,"qp_solver_options" : qpsolver_options, "hessian_approximation": "exact", "UserHM": True}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.input("x0").set([0,1])
      solver.input("lbx").set([-10,-10])
      solver.input("ubx").set([10,10])
      solver.input("lbg").set([0])
      solver.input("ubg").set([1.8])
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],4.64801220074552e-3,6,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],9.318651964592811e-1,5,str(Solver))
      self.assertAlmostEqual(solver.output("x")[1],8.68134821123689e-1,5,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[0],0,8,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[1],0,4,str(Solver))
      self.assertAlmostEqual(solver.output("lam_g")[0],4.75846495145007e-2,5,str(Solver))
      
  def testactiveUBX(self):
    self.message("active UBX")
    x=SX("x")
    y=SX("y")
    
    nlp=SXFunction(nlIn(x=vertcat([x,y])),nlOut(f=(1-x)**2+100*(y-x**2)**2,g=x+y))
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-8,"TolOpti":1e-20,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order","qp_solver": qpsolver,"qp_solver_options" : qpsolver_options, "hessian_approximation": "exact", "UserHM": True}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.input("x0").set([0,1])
      solver.input("lbx").set([-10,0])
      solver.input("ubx").set([10,0.9])
      solver.input("lbg").set([-10])
      solver.input("ubg").set([10])
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],2.626109721583e-3,6,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],9.4882542279172277e-01,6,str(Solver))
      self.assertAlmostEqual(solver.output("x")[1],0.9,6,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[0],0,8,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[1],5.39346608659e-2,4,str(Solver))
      self.assertAlmostEqual(solver.output("lam_g")[0],0,8,str(Solver))
      
  def test_QP(self):
    self.message("QP")

    N = 50 

    x = ssym("x",N)
    x0 = DMatrix(range(N))
    H = diag(range(1,N+1))
    obj = 0.5*mul([(x-x0).T,H,(x-x0)])

    nlp = SXFunction(nlIn(x=x),nlOut(f=obj))
    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-8,"tol_pr":1e-10,"TolOpti":1e-25,"hessian_approximation":"limited-memory","max_iter":100,"maxiter":100, "MaxIter": 100,"print_level":0,"qp_solver": qpsolver,"qp_solver_options" : qpsolver_options}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.input("lbx").setAll(-1000)
      solver.input("ubx").setAll(1000)
      solver.solve()
      self.checkarray(solver.output("x"),x0,str(Solver),digits=2)
      self.assertAlmostEqual(solver.output("f")[0],0,3,str(Solver))
      self.checkarray(solver.output("lam_x"),DMatrix.zeros(N,1),str(Solver),digits=4)
      
      
  def test_tol_pr(self):
    return
    self.message("Low tol_pr")
    H = DMatrix([[1,-1],[-1,2]])
    G = DMatrix([-2,-6])
    A =  DMatrix([[1, 1],[-1, 2],[2, 1]])

    LBA = DMatrix([-inf]*3)
    UBA = DMatrix([2, 2, 3])

    LBX = DMatrix([0.5,0])
    UBX = DMatrix([0.5,inf])

    x=ssym("x",2)
    nlp=SXFunction(nlIn(x=x),nlOut(f=0.5*mul([x.T,H,x])+mul(G.T,x),g=mul(A,x)))

    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-8,"tol_pr":1e-10,"TolOpti":1e-25,"hessian_approximation":"limited-memory","max_iter":100, "maxiter": 100,"MaxIter": 100,"print_level":0,"qp_solver": qpsolver,"qp_solver_options" : qpsolver_options, "fixed_variable_treatment": "make_constraint"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
          
      solver.init()
      solver.input("lbx").set(LBX)
      solver.input("ubx").set(UBX)
      solver.input("lbg").set(LBA)
      solver.input("ubg").set(UBA)

      solver.solve()

      self.assertAlmostEqual(solver.output()[0],0.5,6,str(qpsolver))
      self.assertAlmostEqual(solver.output()[1],1.25,6,str(qpsolver))
    
      self.assertAlmostEqual(solver.output("lam_x")[0],4.75,6,str(qpsolver))
      self.assertAlmostEqual(solver.output("lam_x")[1],0,6,str(qpsolver))

      self.checkarray(solver.output("lam_g"),DMatrix([0,2,0]),str(qpsolver),digits=6)
      
      self.assertAlmostEqual(solver.output("f")[0],-7.4375,6,str(qpsolver))
      
  def test_QP2(self):
    H = DMatrix([[1,-1],[-1,2]])
    G = DMatrix([-2,-6])
    A =  DMatrix([[1, 1],[-1, 2],[2, 1]])

    LBA = DMatrix([-inf]*3)
    UBA = DMatrix([2, 2, 3])

    LBX = DMatrix([0.5,0])
    UBX = DMatrix([0.5,inf])

    x=ssym("x",2)
    nlp=SXFunction(nlIn(x=x),nlOut(f=0.5*mul([x.T,H,x])+mul(G.T,x),g=mul(A,x)))

    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-8,"TolOpti":1e-25,"hessian_approximation":"limited-memory","max_iter":100, "maxiter": 100,"MaxIter": 100,"print_level":0,"qp_solver": qpsolver,"qp_solver_options" : qpsolver_options, "fixed_variable_treatment": "make_constraint"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
          
      solver.init()
      solver.input("lbx").set(LBX)
      solver.input("ubx").set(UBX)
      solver.input("lbg").set(LBA)
      solver.input("ubg").set(UBA)

      solver.solve()

      self.assertAlmostEqual(solver.output()[0],0.5,6,str(qpsolver))
      self.assertAlmostEqual(solver.output()[1],1.25,6,str(qpsolver))
    
      self.assertAlmostEqual(solver.output("lam_x")[0],4.75,6,str(qpsolver))
      self.assertAlmostEqual(solver.output("lam_x")[1],0,6,str(qpsolver))

      self.checkarray(solver.output("lam_g"),DMatrix([0,2,0]),str(qpsolver),digits=6)
      
      self.assertAlmostEqual(solver.output("f")[0],-7.4375,6,str(qpsolver))
      
      solver = Solver(nlp)
      for k,v in ({"tol":1e-8,"TolOpti":1e-25,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "maxiter": 100,"MaxIter": 100,"print_level":0,"qp_solver": qpsolver,"qp_solver_options" : qpsolver_options, "fixed_variable_treatment": "make_constraint"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
          
      solver.init()
      solver.input("lbx").set(LBX)
      solver.input("ubx").set(UBX)
      solver.input("lbg").set(LBA)
      solver.input("ubg").set(UBA)

      solver.solve()

      self.assertAlmostEqual(solver.output()[0],0.5,6,str(qpsolver))
      self.assertAlmostEqual(solver.output()[1],1.25,6,str(qpsolver))
    
      self.assertAlmostEqual(solver.output("lam_x")[0],4.75,6,str(qpsolver))
      self.assertAlmostEqual(solver.output("lam_x")[1],0,6,str(qpsolver))

      self.checkarray(solver.output("lam_g"),DMatrix([0,2,0]),str(qpsolver),digits=6)
      
      self.assertAlmostEqual(solver.output("f")[0],-7.4375,6,str(qpsolver))

  def test_QP2_unconvex(self):
    H = DMatrix([[1,-1],[-1,-2]])
    G = DMatrix([-2,-6])
    A =  DMatrix([[1, 1],[-1, 2],[2, 1]])
    
    LBA = DMatrix([-inf]*3)
    UBA = DMatrix([2, 2, 3])

    LBX = DMatrix([0]*2)
    UBX = DMatrix([inf]*2)

    x=ssym("x",2)
    nlp=SXFunction(nlIn(x=x),nlOut(f=0.5*mul([x.T,H,x])+mul(G.T,x),g=mul(A,x)))

    for Solver in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      for k,v in ({"tol":1e-8,"TolOpti":1e-25,"hessian_approximation":"limited-memory","max_iter":100, "maxiter": 100,"MaxIter": 100,"print_level":0,"qp_solver": qpsolver,"qp_solver_options" : qpsolver_options, "fixed_variable_treatment": "make_constraint"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
          
      solver.init()
      solver.input("lbx").set(LBX)
      solver.input("ubx").set(UBX)
      solver.input("lbg").set(LBA)
      solver.input("ubg").set(UBA)

      solver.solve()

      self.assertAlmostEqual(solver.output()[0],2.0/3,6,str(qpsolver))
      self.assertAlmostEqual(solver.output()[1],4.0/3,6,str(qpsolver))
    
      self.assertAlmostEqual(solver.output("lam_x")[0],0,6,str(qpsolver))
      self.assertAlmostEqual(solver.output("lam_x")[1],0,6,str(qpsolver))

      self.checkarray(solver.output("lam_g"),DMatrix([4+8.0/9,20.0/9,0]),str(qpsolver),digits=6)
      
      self.assertAlmostEqual(solver.output("f")[0],-10-16.0/9,6,str(qpsolver))

      solver = Solver(nlp)
      for k,v in ({"tol":1e-8,"TolOpti":1e-25,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "maxiter": 100,"MaxIter": 100,"print_level":0,"qp_solver": qpsolver,"qp_solver_options" : qpsolver_options, "fixed_variable_treatment": "make_constraint"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
          
      solver.init()
      solver.input("lbx").set(LBX)
      solver.input("ubx").set(UBX)
      solver.input("lbg").set(LBA)
      solver.input("ubg").set(UBA)

      solver.solve()

      self.assertAlmostEqual(solver.output()[0],2.0/3,6,str(qpsolver))
      self.assertAlmostEqual(solver.output()[1],4.0/3,6,str(qpsolver))
    
      self.assertAlmostEqual(solver.output("lam_x")[0],0,6,str(qpsolver))
      self.assertAlmostEqual(solver.output("lam_x")[1],0,6,str(qpsolver))

      self.checkarray(solver.output("lam_g"),DMatrix([4+8.0/9,20.0/9,0]),str(qpsolver),digits=6)
      
      self.assertAlmostEqual(solver.output("f")[0],-10-16.0/9,6,str(qpsolver))
      
if __name__ == '__main__':
    unittest.main()
    print solvers

