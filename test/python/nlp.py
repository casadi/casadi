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
  solvers.append((WorhpSolver,{}))
  print "Will test WorhpSolver"
except:
  pass
  
try:
  solvers.append((IpoptSolver,{}))
  print "Will test IpoptSolver"
except:
  pass

try:
  qpsolver_options = {"nlp_solver": IpoptSolver, "nlp_solver_options": {"tol": 1e-12} }
  solvers.append((SQPMethod,{"qp_solver": NLPQPSolver,"qp_solver_options": qp_solver_options}))
  print "Will test SQPMethod"
except:
  pass

#try:
#  solvers.append(KnitroSolver)
#  print "Will test KnitroSolver"
#except:
#  pass
  

class NLPtests(casadiTestCase):
  def testIPOPT(self):
    x=SX("x")
    nlp=SXFunction(nlpIn(x=x),nlpOut(f=(x-1)**2,g=x))
    
    for Solver, solver_options in solvers:
      self.message("trivial " + str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-5,"hessian_approximation":"limited-memory","max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order" }).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
       
      solver.init()
      solver.setInput([-10],"lbx")
      solver.setInput([10],"ubx")
      solver.setInput([-10],"lbg")
      solver.setInput([10],"ubg")
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,9,str(Solver))
      self.assertAlmostEqual(solver.output("g")[0],1,9,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[0],0,9,str(Solver))
      self.assertAlmostEqual(solver.output("lam_g")[0],0,9,str(Solver))
      
  def testIPOPT_par(self):
    x=SX("x")
    p=SX("p")
    nlp=SXFunction(nlpIn(x=x,p=p),nlpOut(f=(x-p)**2,g=x))
    
    for Solver, solver_options in solvers:
      self.message("trivial " + str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-5,"hessian_approximation":"limited-memory","max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([-10],"lbx")
      solver.setInput([10],"ubx")
      solver.setInput([-10],"lbg")
      solver.setInput([10],"ubg")
      solver.setInput(1,"p")
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,9,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[0],0,9,str(Solver))
      self.assertAlmostEqual(solver.output("lam_g")[0],0,9,str(Solver))
      
  def testIPOPTinf(self):
    self.message("trivial IPOPT, infinity bounds")
    x=SX("x")
    nlp=SXFunction(nlpIn(x=x),nlpOut(f=(x-1)**2,g=x))
    
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-5,"hessian_approximation":"limited-memory","max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([-Inf],"lbx")
      solver.setInput([Inf],"ubx")
      solver.setInput([-Inf],"lbg")
      solver.setInput([Inf],"ubg")

      if 'Worhp' in str(Solver):
        with self.assertRaises(Exception):
          solver.solve()
        return




      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,7,str(Solver) + str(solver.output("x")[0]-1))
      self.assertAlmostEqual(solver.output("lam_x")[0],0,9,str(Solver))
      self.assertAlmostEqual(solver.output("lam_g")[0],0,9,str(Solver))
      
  def testIPOPTrb(self):
    self.message("rosenbrock, limited-memory hessian approx")
    x=SX("x")
    y=SX("y")
    
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=(1-x)**2+100*(y-x**2)**2))
    
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-9,"TolOpti":1e-14,"hessian_approximation":"limited-memory","max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([-10]*2,"lbx")
      solver.setInput([10]*2,"ubx")
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
    
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=(1-x)**2+100*(y-x**2)**2,g=x+y))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-8,"TolOpti":1e-20,"hessian_approximation":"limited-memory","max_iter":1000, "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([-10]*2,"lbx")
      solver.setInput([10]*2,"ubx")
      solver.setInput([-10],"lbg")
      solver.setInput([10],"ubg")
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
    
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=(1-x)**2+100*(y-x**2)**2,g=x+y))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-8,"TolOpti":1e-20,"hessian_approximation":"limited-memory","max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([0,1],"x0")
      solver.setInput([-10,1],"lbx")
      solver.setInput([10,1],"ubx")
      solver.setInput([-10],"lbg")
      solver.setInput([10],"ubg")

      if 'Worhp' in str(Solver):
        with self.assertRaises(Exception):
          solver.solve()
        return




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
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=obj,g=x**2+y**2))
    
    c_r = 4.56748075136258e-02;
    x_r = [7.86415156987791e-01,6.17698316967954e-01]
    
    sigma=SX("sigma")
    lambd=SX("lambd")
    h=SXFunction(hessLagIn(x=vertcat([x,y]),lam_f=sigma,lam_g=lambd),
                 hessLagOut(hess=sigma*hessian(obj,vertcat([x,y]))+lambd*hessian(nlp.outputExpr("g"),vertcat([x,y]))))
    h.init()
    h.setInput([0.5,0.5])
    h.setInput(-40,1)
    h.setInput(1,2)
    h.evaluate()
    print h.output()
    
    solver = None
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      solver.setOption("hess_lag",h)
      for k,v in ({"tol":1e-10,"TolOpti":1e-20,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "MaxIter": 100,"derivative_test":"second-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
          
      solver.init()
      solver.setInput([0.5,0.5],"x0")
      solver.setInput([-10]*2,"lbx")
      solver.setInput([10]*2,"ubx")
      solver.setInput([0],"lbg")
      solver.setInput([1],"ubg")
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
    try:
      solver = IpoptSolver(nlp)
    except:
      return # No IPOPT available
    solver.setOption(solver_options)
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
    solver.setInput([-10]*2,"lbx")
    solver.setInput([10]*2,"ubx")
    solver.setInput([0],"lbg")
    solver.setInput([1],"ubg")
    solver.setInput(oldsolver.output("x"),"x0")
    solver.setInput(oldsolver.output("lam_g"),"lam_g0")
    solver.output("lam_x").set(oldsolver.output("lam_x"))
    
    
    solver.solve()

  def testIPOPTrhb2_gen(self):
    self.message("rosenbrock, exact hessian generated, constrained")
    x=SX("x")
    y=SX("y")
    
    obj = (1-x)**2+100*(y-x**2)**2
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=obj,g=x**2+y**2))
    
    c_r = 4.56748075136258e-02;
    x_r = [7.86415156987791e-01,6.17698316967954e-01]
    
    sigma=SX("sigma")
    lambd=SX("lambd")
  
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-12,"TolOpti":1e-20,"hessian_approximation":"exact","UserHM":True,"max_iter":200, "MaxIter": 100,"print_level":1,"derivative_test":"second-order", "toldx": 1e-15, "tolgl": 1e-15}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
          
      solver.init()
      solver.setInput([0.5,0.5],"x0")
      solver.setInput([-10]*2,"lbx")
      solver.setInput([10]*2,"ubx")
      solver.setInput([0],"lbg")
      solver.setInput([1],"ubg")
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
    nlp=SXFunction(nlpIn(x=vertcat([x,y]),p=p),nlpOut(f=obj,g=x**2+y**2))
    
    c_r = 4.56748075136258e-02;
    x_r = [7.86415156987791e-01,6.17698316967954e-01]
    
    sigma=SX("sigma")
    lambd=SX("lambd")
    h=SXFunction(hessLagIn(x=vertcat([x,y]),lam_f=sigma,lam_g=lambd),
                 hessLagOut(hess=sigma*hessian(obj,vertcat([x,y]))+lambd*hessian(nlp.outputExpr("g"),vertcat([x,y]))))

    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      solver.setOption("hess_lag",h)
      for k,v in ({"tol":1e-10,"TolOpti":1e-20,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "MaxIter": 100,"print_level":1,"derivative_test":"second-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([0.5,0.5],"x0")
      solver.setInput([-10]*2,"lbx")
      solver.setInput([10]*2,"ubx")
      solver.setInput([0],"lbg")
      solver.setInput([1],"ubg")
      solver.setInput([1],"p")
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
    nlp=SXFunction(nlpIn(x=vertcat([x,y]),p=p),nlpOut(f=obj,g=x**2+y**2))
    
    c_r = 4.56748075136258e-02;
    x_r = [7.86415156987791e-01,6.17698316967954e-01]
    
    sigma=SX("sigma")
    lambd=SX("lambd")
  
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-10,"TolOpti":1e-20,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "MaxIter": 100,"print_level":1,"derivative_test":"second-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
          
      solver.init()
      solver.setInput([0.5,0.5],"x0")
      solver.setInput([-10]*2,"lbx")
      solver.setInput([10]*2,"ubx")
      solver.setInput([0],"lbg")
      solver.setInput([1],"ubg")
      solver.setInput([1],"p")
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
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=obj))
    
    sigma=SX("sigma")
    
    h=SXFunction(hessLagIn(x=vertcat([x,y]),lam_f=sigma),
                 hessLagOut(hess=sigma*hessian(obj,vertcat([x,y]))))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      solver.setOption("hess_lag",h)
      for k,v in ({"tol":1e-10,"TolOpti":1e-20,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      #solver.setOption("verbose",True)
      solver.init()
      solver.setInput([-10]*2,"lbx")
      solver.setInput([10]*2,"ubx")
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
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=obj))
    
    sigma=SX("sigma")
    
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-10,"TolOpti":1e-20,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      #solver.setOption("verbose",True)
      solver.init()
      solver.setInput([-10]*2,"lbx")
      solver.setInput([10]*2,"ubx")
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
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=obj))
    
    sigma=SX("sigma")
    
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-10,"TolOpti":1e-20,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      #solver.setOption("verbose",True)
      solver.init()
      solver.setInput([1,-10],"lbx")
      solver.setInput([1,10],"ubx")

      if 'Worhp' in str(Solver):
        with self.assertRaises(Exception):
          solver.solve()
        return



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
    nlp=SXFunction(nlpIn(x=vertcat([x,y]),p=p),nlpOut(f=obj))
    
    sigma=SX("sigma")
    
    h=SXFunction(hessLagIn(x=vertcat([x,y]),p=p,lam_f=sigma),
                 hessLagOut(hess=sigma*hessian(obj,vertcat([x,y]))))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      solver.setOption("hess_lag",h)
      for k,v in ({"tol":1e-10,"TolOpti":1e-20,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      #solver.setOption("verbose",True)
      solver.init()
      solver.setInput([-10]*2,"lbx")
      solver.setInput([10]*2,"ubx")
      solver.setInput(1,"p")
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
    nlp=SXFunction(nlpIn(x=vertcat([x,y]),p=p),nlpOut(f=obj))
    
    sigma=SX("sigma")
    
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-10,"TolOpti":1e-20,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      #solver.setOption("verbose",True)
      solver.init()
      solver.setInput([-10]*2,"lbx")
      solver.setInput([10]*2,"ubx")
      solver.setInput(1,"p")
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
    nlp=MXFunction(nlpIn(x=x),nlpOut(f=norm_2(x-X0),g=2*x))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-8,"max_iter":103, "MaxIter": 103,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([-10]*N,"lbx")
      solver.setInput([10]*N,"ubx")
      solver.setInput([-10]*N,"lbg")
      solver.setInput([10]*N,"ubg")
      solver.solve()
      print "residuals"
      print array(solver.output("x")).squeeze()-x0
      self.assertAlmostEqual(solver.output("f")[0],0,10,str(Solver))
      self.checkarray(array(solver.output("x")).squeeze(),x0,str(Solver),digits=8)
      self.checkarray(solver.output("lam_x"),DMatrix([0]*10),8,str(Solver),digits=8)
      self.assertAlmostEqual(solver.output("lam_g")[1],0,8,str(Solver))
      
  def testIPOPTnoc(self):
    self.message("trivial IPOPT, no constraints")
    """ There is an assertion error thrown, but still it works"""
    x=ssym("x")
    nlp=SXFunction(nlpIn(x=x),nlpOut(f=(x-1)**2))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-10,"max_iter":103, "MaxIter": 103,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver = IpoptSolver(nlp)
      solver.init()
      solver.setInput([-10],"lbx")
      solver.setInput([10],"ubx")
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,9,str(Solver))
    
  def testIPOPTmx(self):
    self.message("trivial IPOPT, using MX")
    x=MX("x")
    nlp=MXFunction(nlpIn(x=x),nlpOut(f=(x-1)**2,g=2*x))
    
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-10,"max_iter":103, "MaxIter": 103,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([-10],"lbx")
      solver.setInput([10],"ubx")
      solver.setInput([-10],"lbg")
      solver.setInput([10],"ubg")
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,9,str(Solver))
    
  def testIPOPTc(self):
    self.message("trivial, overconstrained")
    x=SX("x")
    nlp=SXFunction(nlpIn(x=x),nlpOut(f=(x-1)**2,g=vertcat([x,x,x])))
    
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-5,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([-10],"lbx")
      solver.setInput([10],"ubx")
      solver.setInput([-10, -10, -10],"lbg")
      solver.setInput([10, 10, 10],"ubg")
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],0,9,str(Solver) )
      self.assertAlmostEqual(solver.output("x")[0],1,5,str(Solver))
    
  def testIPOPTc2(self):
    self.message("trivial2, overconstrained")
    x=SX("x")
    nlp=SXFunction(nlpIn(x=x),nlpOut(f=(x-1)**2,g=vertcat([x,x,x+x])))
    
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-10,"max_iter":100, "hessian_approximation": "limited-memory", "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([-10],"lbx")
      solver.setInput([10],"ubx")
      solver.setInput([-10, -10, -10],"lbg")
      solver.setInput([10, 10, 10],"ubg")
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,8,str(Solver))
    
  def testIPOPTcmx(self):
    self.message("trivial , overconstrained, using MX")
    x=MX("x")
    nlp=MXFunction(nlpIn(x=x),nlpOut(f=(x-1)**2,g=vertcat([2*x,3*x,4*x])))
    
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-10,"max_iter":100, "hessian_approximation": "limited-memory", "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([-10],"lbx")
      solver.setInput([10],"ubx")
      solver.setInput([-10,-10,-10],"lbg")
      solver.setInput([10,10,10],"ubg")
      solver.solve()
      self.assertAlmostEqual(solver.output("f")[0],0,9,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,8,str(Solver))

  def testIPOPTdeg(self):
    self.message("degenerate optimization IPOPT")
    x=SX("x")
    y=SX("y")
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=0,g=vertcat([x-y,x])))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-5,"max_iter":100, "hessian_approximation": "limited-memory", "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([-10, -10],"lbx")
      solver.setInput([10, 10],"ubx")
      solver.setInput([0, 3],"lbg")
      solver.setInput([0, 3],"ubg")
      solver.solve()
      self.assertAlmostEqual(solver.output("x")[0],solver.output("x")[1],10,"IPOPT")

  def testIPOPTdegc(self):
    self.message("degenerate optimization IPOPT, overconstrained")
    x=SX("x")
    y=SX("y")
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=0,g=vertcat([x-y,x,x+y])))
    
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-5,"max_iter":100, "hessian_approximation": "limited-memory", "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)

      solver.init()
      solver.setInput([-10, -10],"lbx")
      solver.setInput([10, 10],"ubx")
      solver.setInput([0, 3 , -10],"lbg")
      solver.setInput([0, 3, 10],"ubg")
      solver.solve()
      # todo: catch error when set([0, 3 , 5]) two times
      self.assertAlmostEqual(solver.output("x")[0],solver.output("x")[1],10,"IPOPT")
      
  def testXfreeChange(self):
    self.message("Change in X settings")
    x=SX("x")
    y=SX("y")
    
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=(1-x)**2+100*(y-x**2)**2,g=x+y))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-8,"TolOpti":1e-20,"hessian_approximation":"limited-memory","max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([0,1],"x0")
      solver.setInput([-10,-10],"lbx")
      solver.setInput([10,10],"ubx")
      solver.setInput([-10],"lbg")
      solver.setInput([10],"ubg")
      solver.solve()
      solver.setInput([-10,1],"lbx")
      solver.setInput([10,1],"ubx")
      solver.setInput([-10],"lbg")
      solver.setInput([10],"ubg")

      if 'Worhp' in str(Solver):
        with self.assertRaises(Exception):
          solver.solve()
        return


      solver.solve()
      
      self.assertAlmostEqual(solver.output("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.output("x")[0],1,7,str(Solver))
      self.assertAlmostEqual(solver.output("x")[1],1,7,str(Solver))

  def testactiveLBX(self):
    self.message("active LBX")
    x=SX("x")
    y=SX("y")
    
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=(1-x)**2+100*(y-x**2)**2,g=x+y))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-8,"TolOpti":1e-20,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order", "hessian_approximation": "exact", "UserHM": True}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([0,1],"x0")
      solver.setInput([-10,1.2],"lbx")
      solver.setInput([10,2],"ubx")
      solver.setInput([-10],"lbg")
      solver.setInput([10],"ubg")
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
    
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=(1-x)**2+100*(y-x**2)**2,g=x+y))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-8,"TolOpti":1e-20,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order", "hessian_approximation": "exact", "UserHM": True }).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([0,1],"x0")
      solver.setInput([-10,-10],"lbx")
      solver.setInput([10,10],"ubx")
      solver.setInput([2.2],"lbg")
      solver.setInput([10],"ubg")
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
    
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=(1-x)**2+100*(y-x**2)**2,g=x+y))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-8,"TolOpti":1e-20,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order", "hessian_approximation": "exact", "UserHM": True}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([0,1],"x0")
      solver.setInput([-10,-10],"lbx")
      solver.setInput([10,10],"ubx")
      solver.setInput([0],"lbg")
      solver.setInput([1.8],"ubg")
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
    
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=(1-x)**2+100*(y-x**2)**2,g=x+y))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-8,"TolOpti":1e-20,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order", "hessian_approximation": "exact", "UserHM": True}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([0,1],"x0")
      solver.setInput([-10,0],"lbx")
      solver.setInput([10,0.9],"ubx")
      solver.setInput([-10],"lbg")
      solver.setInput([10],"ubg")
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

    nlp = SXFunction(nlpIn(x=x),nlpOut(f=obj))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-8,"tol_pr":1e-10,"TolOpti":1e-25,"hessian_approximation":"limited-memory","max_iter":100, "MaxIter": 100,"print_level":0}).iteritems():
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
    nlp=SXFunction(nlpIn(x=x),nlpOut(f=0.5*mul([x.T,H,x])+mul(G.T,x),g=mul(A,x)))

    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-8,"tol_pr":1e-10,"TolOpti":1e-25,"hessian_approximation":"limited-memory","max_iter":100,"MaxIter": 100,"print_level":0, "fixed_variable_treatment": "make_constraint"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
          
      solver.init()
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lbg")
      solver.setInput(UBA,"ubg")

      solver.solve()

      self.assertAlmostEqual(solver.output()[0],0.5,6,str(Solver))
      self.assertAlmostEqual(solver.output()[1],1.25,6,str(Solver))
    
      self.assertAlmostEqual(solver.output("lam_x")[0],4.75,6,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[1],0,6,str(Solver))

      self.checkarray(solver.output("lam_g"),DMatrix([0,2,0]),str(Solver),digits=6)
      
      self.assertAlmostEqual(solver.output("f")[0],-7.4375,6,str(Solver))
      
  def test_QP2(self):
    H = DMatrix([[1,-1],[-1,2]])
    G = DMatrix([-2,-6])
    A =  DMatrix([[1, 1],[-1, 2],[2, 1]])

    LBA = DMatrix([-inf]*3)
    UBA = DMatrix([2, 2, 3])

    LBX = DMatrix([0.5,0])
    UBX = DMatrix([0.5,inf])

    x=ssym("x",2)
    nlp=SXFunction(nlpIn(x=x),nlpOut(f=0.5*mul([x.T,H,x])+mul(G.T,x),g=mul(A,x)))

    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-8,"TolOpti":1e-25,"hessian_approximation":"limited-memory","max_iter":100,"MaxIter": 100,"print_level":0, "fixed_variable_treatment": "make_constraint"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
          
      solver.init()
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lbg")
      solver.setInput(UBA,"ubg")
      if 'Worhp' in str(Solver):
        with self.assertRaises(Exception):
          solver.solve()
        return

      solver.solve()

      self.assertAlmostEqual(solver.output()[0],0.5,6,str(Solver))
      self.assertAlmostEqual(solver.output()[1],1.25,6,str(Solver))
    
      self.assertAlmostEqual(solver.output("lam_x")[0],4.75,6,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[1],0,6,str(Solver))

      self.checkarray(solver.output("lam_g"),DMatrix([0,2,0]),str(Solver),digits=6)
      
      self.assertAlmostEqual(solver.output("f")[0],-7.4375,6,str(Solver))
      
      solver = Solver(nlp)
      for k,v in ({"tol":1e-8,"TolOpti":1e-25,"hessian_approximation":"exact","UserHM":True,"max_iter":100,"MaxIter": 100,"print_level":0, "fixed_variable_treatment": "make_constraint"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
          
      solver.init()
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lbg")
      solver.setInput(UBA,"ubg")

      solver.solve()

      self.assertAlmostEqual(solver.output()[0],0.5,6,str(Solver))
      self.assertAlmostEqual(solver.output()[1],1.25,6,str(Solver))
    
      self.assertAlmostEqual(solver.output("lam_x")[0],4.75,6,str(Solver))
      self.assertAlmostEqual(solver.output("lam_x")[1],0,6,str(Solver))

      self.checkarray(solver.output("lam_g"),DMatrix([0,2,0]),str(Solver),digits=6)
      
      self.assertAlmostEqual(solver.output("f")[0],-7.4375,6,str(Solver))

  def test_QP2_unconvex(self):
    H = DMatrix([[1,-1],[-1,-2]])
    G = DMatrix([-2,-6])
    A =  DMatrix([[1, 1],[-1, 2],[2, 1]])
    
    LBA = DMatrix([-inf]*3)
    UBA = DMatrix([2, 2, 3])

    LBX = DMatrix([0]*2)
    UBX = DMatrix([inf]*2)

    x=ssym("x",2)
    nlp=SXFunction(nlpIn(x=x),nlpOut(f=0.5*mul([x.T,H,x])+mul(G.T,x),g=mul(A,x)))

    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-8,"TolOpti":1e-25,"hessian_approximation":"limited-memory","max_iter":100,"MaxIter": 100,"print_level":0, "fixed_variable_treatment": "make_constraint"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
          
      solver.init()
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lbg")
      solver.setInput(UBA,"ubg")

      solver.solve()

      self.assertAlmostEqual(solver.output()[0],2.0/3,6,str(solver))
      self.assertAlmostEqual(solver.output()[1],4.0/3,6,str(solver))
    
      self.assertAlmostEqual(solver.output("lam_x")[0],0,6,str(solver))
      self.assertAlmostEqual(solver.output("lam_x")[1],0,6,str(solver))

      self.checkarray(solver.output("lam_g"),DMatrix([4+8.0/9,20.0/9,0]),str(solver),digits=6)
      
      self.assertAlmostEqual(solver.output("f")[0],-10-16.0/9,6,str(solver))

      solver = Solver(nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-8,"TolOpti":1e-25,"hessian_approximation":"exact","UserHM":True,"max_iter":100,"MaxIter": 100,"print_level":0, "fixed_variable_treatment": "make_constraint"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
          
      solver.init()
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lbg")
      solver.setInput(UBA,"ubg")

      solver.solve()

      self.assertAlmostEqual(solver.output()[0],2.0/3,6,str(solver))
      self.assertAlmostEqual(solver.output()[1],4.0/3,6,str(solver))
    
      self.assertAlmostEqual(solver.output("lam_x")[0],0,6,str(solver))
      self.assertAlmostEqual(solver.output("lam_x")[1],0,6,str(solver))

      self.checkarray(solver.output("lam_g"),DMatrix([4+8.0/9,20.0/9,0]),str(solver),digits=6)
      
      self.assertAlmostEqual(solver.output("f")[0],-10-16.0/9,6,str(solver))
      
if __name__ == '__main__':
    unittest.main()
    print solvers

