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
from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *
import itertools

#CasadiOptions.setCatchErrorsPython(False)

solvers= []
 
try:
  NlpSolver.loadPlugin("worhp")
  solvers.append(("worhp",{}))
  print "Will test worhp"
except:
  pass
  
try:
  NlpSolver.loadPlugin("ipopt")
  solvers.append(("ipopt",{}))
  print "Will test ipopt"
except:
  pass

try:
  NlpSolver.loadPlugin("snopt")
  solvers.append(("snopt",{"Verify level": 3,"detect_linear": True,"Major optimality tolerance":1e-12,"Minor feasibility tolerance":1e-12,"Major feasibility tolerance":1e-12}))
  print "Will test snopt"
except:
  pass

try:
  NlpSolver.loadPlugin("ipopt")
  NlpSolver.loadPlugin("sqpmethod")
  qp_solver_options = {"nlp_solver": "ipopt", "nlp_solver_options": {"tol": 1e-12} }
  solvers.append(("sqpmethod",{"qp_solver": "nlp","qp_solver_options": qp_solver_options}))
  print "Will test sqpmethod"
except:
  pass
  
try:
  NlpSolver.loadPlugin("ipopt")
  NlpSolver.loadPlugin("stabilizedsqp")
  qp_solver_options = {"nlp_solver": "ipopt", "nlp_solver_options": {"tol": 1e-12, "print_level": 0, "print_time": False} }
  solvers.append(("stabilizedsqp",{"tol_pr": 1e-9, "tol_du": 1e-9,"stabilized_qp_solver": "qp", "stabilized_qp_solver_options": {"qp_solver": "nlp", "qp_solver_options": qp_solver_options}}))
  print "Will test stabilizedsqp"
except:
  pass
  
try:
  qp_solver_options = {}
  QpSolver.loadPlugin("sqic")
  NlpSolver.loadPlugin("stabilizedsqp")
  solvers.append(("stabilizedsqp",{"tol_pr": 1e-9, "tol_du": 1e-9,"stabilized_qp_solver": "qp", "stabilized_qp_solver_options": {"qp_solver": "sqic"}}))
  print "Will test stabilizedsqp"
except:
  pass

"""
try:
  NlpSolver.loadPlugin("knitro")
  solvers.append(("knitro",{}))
  print "Will test knitro"
except:
  pass
"""

class NLPtests(casadiTestCase):

  def testboundsviol(self):
    x=SX.sym("x")
    nlp=SXFunction(nlpIn(x=x),nlpOut(f=(x-1)**2,g=x))
    
    for Solver, solver_options in solvers:
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-5,"hessian_approximation":"limited-memory","max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order" }).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
       
      solver.init()
      solver.setInput([-10],"lbx")
      solver.setInput([-20],"ubx")
      solver.setInput([-10],"lbg")
      solver.setInput([10],"ubg")
      with self.assertRaises(Exception):
        solver.evaluate()

    for Solver, solver_options in solvers:
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-5,"hessian_approximation":"limited-memory","max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order" }).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
       
      solver.init()
      solver.setInput([-10],"lbx")
      solver.setInput([10],"ubx")
      solver.setInput([-10],"lbg")
      solver.setInput([-20],"ubg")
      with self.assertRaises(Exception):
        solver.evaluate()
        
  def testIPOPT(self):
    x=SX.sym("x")
    nlp=SXFunction(nlpIn(x=x),nlpOut(f=(x-1)**2,g=x))
    
    for Solver, solver_options in solvers:
      self.message("trivial " + str(Solver))
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-5,"hessian_approximation":"limited-memory","max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order" }).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
       
      solver.init()
      solver.setInput([-10],"lbx")
      solver.setInput([10],"ubx")
      solver.setInput([-10],"lbg")
      solver.setInput([10],"ubg")
      solver.evaluate()
      self.assertAlmostEqual(solver.getOutput("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[0],1,9,str(Solver))
      self.assertAlmostEqual(solver.getOutput("g")[0],1,9,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,9,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_g")[0],0,9,str(Solver))
      
  def testIPOPT_par(self):
    x=SX.sym("x")
    p=SX.sym("p")
    nlp=SXFunction(nlpIn(x=x,p=p),nlpOut(f=(x-p)**2,g=x))
    
    for Solver, solver_options in solvers:
      self.message("trivial " + str(Solver))
      solver = NlpSolver(Solver, nlp)
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
      solver.evaluate()
      self.assertAlmostEqual(solver.getOutput("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[0],1,9,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,9,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_g")[0],0,9,str(Solver))
      
  def testIPOPTinf(self):
    self.message("trivial IPOPT, infinity bounds")
    x=SX.sym("x")
    nlp=SXFunction(nlpIn(x=x),nlpOut(f=(x-1)**2,g=x))
    
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-5,"hessian_approximation":"limited-memory","max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([-Inf],"lbx")
      solver.setInput([Inf],"ubx")
      solver.setInput([-Inf],"lbg")
      solver.setInput([Inf],"ubg")

      if Solver in ("worhp","knitro"):
        with self.assertRaises(Exception):
          solver.evaluate()
        return




      solver.evaluate()
      self.assertAlmostEqual(solver.getOutput("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[0],1,7,str(Solver) + str(solver.getOutput("x")[0]-1))
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,9,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_g")[0],0,9,str(Solver))
      
  def testIPOPTrb(self):
    self.message("rosenbrock, limited-memory hessian approx")
    x=SX.sym("x")
    y=SX.sym("y")
    
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=(1-x)**2+100*(y-x**2)**2))
    
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-9,"TolOpti":1e-14,"hessian_approximation":"limited-memory","max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([-10]*2,"lbx")
      solver.setInput([10]*2,"ubx")
      solver.evaluate()
      self.assertAlmostEqual(solver.getOutput("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[0],1,6,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[1],1,6,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,5,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,5,str(Solver))
    
  def testIPOPTrb2(self):
    self.message("rosenbrock, limited-memory hessian approx")
    x=SX.sym("x")
    y=SX.sym("y")
    
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=(1-x)**2+100*(y-x**2)**2,g=x+y))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-8,"TolOpti":1e-20,"hessian_approximation":"limited-memory","max_iter":1000, "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([-10]*2,"lbx")
      solver.setInput([10]*2,"ubx")
      solver.setInput([-10],"lbg")
      solver.setInput([10],"ubg")
      solver.evaluate()
      
      digits = 6

      self.assertAlmostEqual(solver.getOutput("f")[0],0,digits,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[0],1,digits,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[1],1,digits,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,5,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,5,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_g")[0],0,5,str(Solver))
      
  def testIPOPTrbf(self):
    self.message("rosenbrock fixed, limited-memory hessian approx")
    x=SX.sym("x")
    y=SX.sym("y")
    
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=(1-x)**2+100*(y-x**2)**2,g=x+y))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
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

      if 'worhp' in str(Solver):
        with self.assertRaises(Exception):
          solver.evaluate()
        return




      solver.evaluate()
      self.assertAlmostEqual(solver.getOutput("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[0],1,7,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[1],1,7,str(Solver))
      if "stabilizedsqp" not in str(Solver):
        self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,6,str(Solver))
        self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,6,str(Solver))
        self.assertAlmostEqual(solver.getOutput("lam_g")[0],0,6,str(Solver))
      
  def test_IPOPTrhb2(self):
    self.message("rosenbrock, exact hessian, constrained")
    x=SX.sym("x")
    y=SX.sym("y")
    
    obj = (1-x)**2+100*(y-x**2)**2
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=obj,g=x**2+y**2))
    
    c_r = 4.56748075136258e-02;
    x_r = [7.86415156987791e-01,6.17698316967954e-01]
    
    sigma=SX.sym("sigma")
    lambd=SX.sym("lambd")
    h=SXFunction(hessLagIn(x=vertcat([x,y]),lam_f=sigma,lam_g=lambd),
                 hessLagOut(hess=sigma*hessian(obj,vertcat([x,y]))+lambd*hessian(nlp.outputExpr("g"),vertcat([x,y]))))
    h.init()
    h.setInput([0.5,0.5])
    h.setInput(-40,1)
    h.setInput(1,2)
    h.evaluate()
    print h.getOutput()
    
    solver = None
    for Solver, solver_options in solvers:
      self.message(Solver)
      solver = NlpSolver(Solver, nlp)
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
      solver.evaluate()
      
      digits = 5
        
      self.assertAlmostEqual(solver.getOutput("f")[0],c_r,digits,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[0],x_r[0],digits,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[1],x_r[1],digits,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,8,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,8,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_g")[0],0.12149655447670,6,str(Solver))
      
  def test_warmstart(self):
  
    x=SX.sym("x")
    y=SX.sym("y")
    
    obj = (1-x)**2+100*(y-x**2)**2
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=obj,g=x**2+y**2))
    
    c_r = 4.56748075136258e-02;
    x_r = [7.86415156987791e-01,6.17698316967954e-01]
    
    for Solver, solver_options in solvers:
      self.message(Solver)
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-10,"TolOpti":1e-20,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "MaxIter": 100,"derivative_test":"second-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([0.5,0.5],"x0")
      solver.setInput([-10]*2,"lbx")
      solver.setInput([10]*2,"ubx")
      solver.setInput([0],"lbg")
      solver.setInput([1],"ubg")
      solver.evaluate()
      
      digits = 5
        
      self.assertAlmostEqual(solver.getOutput("f")[0],c_r,digits,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[0],x_r[0],digits,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[1],x_r[1],digits,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,8,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,8,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_g")[0],0.12149655447670,6,str(Solver))

      self.message(":warmstart")
      if "ipopt" in str(Solver):
        oldsolver=solver
        solver = NlpSolver(Solver, nlp)
        solver.setOption(solver_options)
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
        solver.setInput(oldsolver.getOutput("x"),"x0")
        solver.setInput(oldsolver.getOutput("lam_g"),"lam_g0")
        solver.setOutput(oldsolver.getOutput("lam_x"),"lam_x")
        
        
        solver.evaluate()

  def testIPOPTrhb2_gen(self):
    self.message("rosenbrock, exact hessian generated, constrained")
    x=SX.sym("x")
    y=SX.sym("y")
    
    obj = (1-x)**2+100*(y-x**2)**2
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=obj,g=x**2+y**2))
    
    c_r = 4.56748075136258e-02;
    x_r = [7.86415156987791e-01,6.17698316967954e-01]
    
    sigma=SX.sym("sigma")
    lambd=SX.sym("lambd")
  
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
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
      solver.evaluate()
      
      digits = 5
      
      self.assertAlmostEqual(solver.getOutput("f")[0],c_r,digits,str(Solver) + str(solver.getOutput("f")[0]) + ":" + str(c_r))
      self.assertAlmostEqual(solver.getOutput("x")[0],x_r[0],digits,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[1],x_r[1],digits,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,8,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,8,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_g")[0],0.12149655447670,6,str(Solver))
      
      
  def test_jacG_empty(self):
    x=SX.sym("x")
    y=SX.sym("y")
    
    obj = (1-x)**2+100*(y-x**2)**2
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=obj,g=1))
    
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      if "worhp"==Solver:
        continue
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      solver.init()
      solver.setInput([0.5,0.5],"x0")
      solver.setInput([-10]*2,"lbx")
      solver.setInput([10]*2,"ubx")
      solver.setInput([0],"lbg")
      solver.setInput([2],"ubg")
      solver.evaluate()
      
      digits = 5
        
      self.checkarray(solver.getOutput("f"),DMatrix([0]),str(Solver),digits=digits)
      self.checkarray(solver.getOutput("x"),DMatrix([1,1]),str(Solver),digits=digits)
      self.checkarray(solver.getOutput("lam_x"),DMatrix([0,0]),str(Solver),digits=digits)
      self.checkarray(solver.getOutput("lam_g"),DMatrix([0]),str(Solver),digits=digits)

  def testIPOPTrhb2_par(self):
    self.message("rosenbrock, exact hessian, constrained, ")
    x=SX.sym("x")
    y=SX.sym("y")
    p=SX.sym("p")
    
    obj = (p-x)**2+100*(y-x**2)**2
    nlp=SXFunction(nlpIn(x=vertcat([x,y]),p=p),nlpOut(f=obj,g=x**2+y**2))
    
    c_r = 4.56748075136258e-02;
    x_r = [7.86415156987791e-01,6.17698316967954e-01]
    
    sigma=SX.sym("sigma")
    lambd=SX.sym("lambd")
    h=SXFunction(hessLagIn(x=vertcat([x,y]),lam_f=sigma,lam_g=lambd,p=p),
                 hessLagOut(hess=sigma*hessian(obj,vertcat([x,y]))+lambd*hessian(nlp.outputExpr("g"),vertcat([x,y]))))

    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
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
      solver.evaluate()
      
      digits = 5
        
      self.assertAlmostEqual(solver.getOutput("f")[0],c_r,digits,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[0],x_r[0],digits,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[1],x_r[1],digits,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,8,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,8,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_g")[0],0.12149655447670,6,str(Solver))

  def testIPOPTrhb2_gen_par(self):
    self.message("rosenbrock, exact hessian generated, constrained, parametric")
    x=SX.sym("x")
    y=SX.sym("y")
    p=SX.sym("p")
    
    obj = (p-x)**2+100*(y-x**2)**2
    nlp=SXFunction(nlpIn(x=vertcat([x,y]),p=p),nlpOut(f=obj,g=x**2+y**2))
    
    c_r = 4.56748075136258e-02;
    x_r = [7.86415156987791e-01,6.17698316967954e-01]
    
    sigma=SX.sym("sigma")
    lambd=SX.sym("lambd")
  
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
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
      solver.evaluate()
      
      digits = 5

      self.assertAlmostEqual(solver.getOutput("f")[0],c_r,digits,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[0],x_r[0],digits,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[1],x_r[1],digits,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,8,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,8,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_g")[0],0.12149655447670,6,str(Solver))
      
  def test_IPOPTrhb(self):
    self.message("rosenbrock, exact hessian")
    x=SX.sym("x")
    y=SX.sym("y")
    
    obj=(1-x)**2+100*(y-x**2)**2
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=obj))
    
    sigma=SX.sym("sigma")
    
    h=SXFunction(hessLagIn(x=vertcat([x,y]),lam_f=sigma),
                 hessLagOut(hess=sigma*hessian(obj,vertcat([x,y]))))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      solver.setOption("hess_lag",h)
      for k,v in ({"tol":1e-10,"TolOpti":1e-20,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      #solver.setOption("verbose",True)
      solver.init()
      solver.setInput([-10]*2,"lbx")
      solver.setInput([10]*2,"ubx")
      solver.evaluate()
      self.assertAlmostEqual(solver.getOutput("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[0],1,9,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[1],1,9,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,8,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,8,str(Solver))

  def testIPOPTrhb_gen(self):
    self.message("rosenbrock, exact hessian generated")
    x=SX.sym("x")
    y=SX.sym("y")
    
    obj=(1-x)**2+100*(y-x**2)**2
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=obj))
    
    sigma=SX.sym("sigma")
    
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-10,"TolOpti":1e-20,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      #solver.setOption("verbose",True)
      solver.init()
      solver.setInput([-10]*2,"lbx")
      solver.setInput([10]*2,"ubx")
      solver.evaluate()
      self.assertAlmostEqual(solver.getOutput("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[0],1,9,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[1],1,9,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,8,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,8,str(Solver))

  def testIPOPTrhb_gen_xnonfree(self):
    self.message("rosenbrock, exact hessian generated, non-free x")
    x=SX.sym("x")
    y=SX.sym("y")
    
    obj=(1-x)**2+100*(y-x**2)**2
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=obj))
    
    sigma=SX.sym("sigma")
    
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-10,"TolOpti":1e-20,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      #solver.setOption("verbose",True)
      solver.init()
      solver.setInput([1,-10],"lbx")
      solver.setInput([1,10],"ubx")

      if 'worhp' in str(Solver):
        with self.assertRaises(Exception):
          solver.evaluate()
        return



      solver.evaluate()
      self.assertAlmostEqual(solver.getOutput("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[0],1,9,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[1],1,6,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,6,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,6,str(Solver))
      
  def testIPOPTrhb_par(self):
    self.message("rosenbrock, exact hessian, parametric")
    x=SX.sym("x")
    y=SX.sym("y")
    
    p=SX.sym("p")
    obj=(p-x)**2+100*(y-x**2)**2
    nlp=SXFunction(nlpIn(x=vertcat([x,y]),p=p),nlpOut(f=obj))
    
    sigma=SX.sym("sigma")
    
    h=SXFunction(hessLagIn(x=vertcat([x,y]),p=p,lam_f=sigma),
                 hessLagOut(hess=sigma*hessian(obj,vertcat([x,y]))))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
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
      solver.evaluate()
      self.assertAlmostEqual(solver.getOutput("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[0],1,9,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[1],1,9,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,8,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,8,str(Solver))

  def testIPOPTrhb_gen_par(self):
    self.message("rosenbrock, exact hessian generated, parametric")
    x=SX.sym("x")
    y=SX.sym("y")
    
    p=SX.sym("p")
    obj=(p-x)**2+100*(y-x**2)**2
    nlp=SXFunction(nlpIn(x=vertcat([x,y]),p=p),nlpOut(f=obj))
    
    sigma=SX.sym("sigma")
    
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-10,"TolOpti":1e-20,"hessian_approximation":"exact","UserHM":True,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      #solver.setOption("verbose",True)
      solver.init()
      solver.setInput([-10]*2,"lbx")
      solver.setInput([10]*2,"ubx")
      solver.setInput(1,"p")
      solver.evaluate()
      self.assertAlmostEqual(solver.getOutput("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[0],1,9,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[1],1,9,str(Solver))
      
  def testIPOPTnorm(self):
    self.message("IPOPT min ||x||^2_2")
    def norm_2(mx):
      return inner_prod(mx,mx)
    N=10
    x=MX.sym("x",N)
    x0=linspace(0,1,N)
    X0=MX(x0)
    nlp=MXFunction(nlpIn(x=x),nlpOut(f=norm_2(x-X0),g=2*x))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-8,"max_iter":103, "MaxIter": 103,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([-10]*N,"lbx")
      solver.setInput([10]*N,"ubx")
      solver.setInput([-10]*N,"lbg")
      solver.setInput([10]*N,"ubg")
      solver.evaluate()
      print "residuals"
      print array(solver.getOutput("x")).squeeze()-x0
      print "bazmeg", solver.getOutput("f")
      self.assertAlmostEqual(solver.getOutput("f")[0],0,10,str(Solver))
      self.checkarray(array(solver.getOutput("x")).squeeze(),x0,str(Solver),digits=8)
      self.checkarray(solver.getOutput("lam_x"),DMatrix([0]*10),8,str(Solver),digits=8)
      self.assertAlmostEqual(solver.getOutput("lam_g")[1],0,8,str(Solver))
      
  def testIPOPTnoc(self):
    self.message("trivial IPOPT, no constraints")
    """ There is an assertion error thrown, but still it works"""
    x=SX.sym("x")
    nlp=SXFunction(nlpIn(x=x),nlpOut(f=(x-1)**2))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-10,"max_iter":103, "MaxIter": 103,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver = NlpSolver("ipopt", nlp)
      solver.init()
      solver.setInput([-10],"lbx")
      solver.setInput([10],"ubx")
      solver.evaluate()
      self.assertAlmostEqual(solver.getOutput("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[0],1,9,str(Solver))
    
  def testIPOPTmx(self):
    self.message("trivial IPOPT, using MX")
    x=MX.sym("x")
    nlp=MXFunction(nlpIn(x=x),nlpOut(f=(x-1)**2,g=2*x))
    
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-10,"max_iter":103, "MaxIter": 103,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([-10],"lbx")
      solver.setInput([10],"ubx")
      solver.setInput([-10],"lbg")
      solver.setInput([10],"ubg")
      solver.evaluate()
      self.assertAlmostEqual(solver.getOutput("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[0],1,9,str(Solver))
    
  def testIPOPTc(self):
    self.message("trivial, overconstrained")
    x=SX.sym("x")
    nlp=SXFunction(nlpIn(x=x),nlpOut(f=(x-1)**2,g=vertcat([x,x,x])))
    
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-5,"max_iter":100, "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([-10],"lbx")
      solver.setInput([10],"ubx")
      solver.setInput([-10, -10, -10],"lbg")
      solver.setInput([10, 10, 10],"ubg")
      solver.evaluate()
      self.assertAlmostEqual(solver.getOutput("f")[0],0,9,str(Solver) )
      self.assertAlmostEqual(solver.getOutput("x")[0],1,5,str(Solver))
    
  def testIPOPTc2(self):
    self.message("trivial2, overconstrained")
    x=SX.sym("x")
    nlp=SXFunction(nlpIn(x=x),nlpOut(f=(x-1)**2,g=vertcat([x,x,x+x])))
    
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-10,"max_iter":100, "hessian_approximation": "limited-memory", "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([-10],"lbx")
      solver.setInput([10],"ubx")
      solver.setInput([-10, -10, -10],"lbg")
      solver.setInput([10, 10, 10],"ubg")
      solver.evaluate()
      self.assertAlmostEqual(solver.getOutput("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[0],1,8,str(Solver))
    
  def testIPOPTcmx(self):
    self.message("trivial , overconstrained, using MX")
    x=MX.sym("x")
    nlp=MXFunction(nlpIn(x=x),nlpOut(f=(x-1)**2,g=vertcat([2*x,3*x,4*x])))
    
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-10,"max_iter":100, "hessian_approximation": "limited-memory", "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([-10],"lbx")
      solver.setInput([10],"ubx")
      solver.setInput([-10,-10,-10],"lbg")
      solver.setInput([10,10,10],"ubg")
      solver.evaluate()
      self.assertAlmostEqual(solver.getOutput("f")[0],0,9,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[0],1,8,str(Solver))

  def testIPOPTdeg(self):
    self.message("degenerate optimization IPOPT")
    x=SX.sym("x")
    y=SX.sym("y")
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=0,g=vertcat([x-y,x])))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-5,"max_iter":100, "hessian_approximation": "limited-memory", "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([-10, -10],"lbx")
      solver.setInput([10, 10],"ubx")
      solver.setInput([0, 3],"lbg")
      solver.setInput([0, 3],"ubg")
      solver.evaluate()
      self.assertAlmostEqual(solver.getOutput("x")[0],solver.getOutput("x")[1],4 if "sqic" in str(solver_options) else 10,"IPOPT")

  def testIPOPTdegc(self):
    self.message("degenerate optimization IPOPT, overconstrained")
    x=SX.sym("x")
    y=SX.sym("y")
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=0,g=vertcat([x-y,x,x+y])))
    
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-5,"max_iter":100, "hessian_approximation": "limited-memory", "MaxIter": 100,"print_level":0,"derivative_test":"first-order"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)

      solver.init()
      solver.setInput([-10, -10],"lbx")
      solver.setInput([10, 10],"ubx")
      solver.setInput([0, 3 , -10],"lbg")
      solver.setInput([0, 3, 10],"ubg")
      solver.evaluate()
      # todo: catch error when set([0, 3 , 5]) two times
      self.assertAlmostEqual(solver.getOutput("x")[0],solver.getOutput("x")[1],4 if "sqic" in str(solver_options) else 10,"IPOPT")
      
  def testXfreeChange(self):
    self.message("Change in X settings")
    x=SX.sym("x")
    y=SX.sym("y")
    
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=(1-x)**2+100*(y-x**2)**2,g=x+y))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
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
      solver.evaluate()
      solver.setInput([-10,1],"lbx")
      solver.setInput([10,1],"ubx")
      solver.setInput([-10],"lbg")
      solver.setInput([10],"ubg")

      if 'worhp' in str(Solver):
        with self.assertRaises(Exception):
          solver.evaluate()
        return


      solver.evaluate()
      
      self.assertAlmostEqual(solver.getOutput("f")[0],0,10,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[0],1,7,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[1],1,7,str(Solver))

  def test_activeLBX(self):
    self.message("active LBX")
    x=SX.sym("x")
    y=SX.sym("y")
    
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=(1-x)**2+100*(y-x**2)**2,g=x+y))
    for Solver, solver_options in solvers:
      self.message(Solver)
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-8,"TolOpti":1e-20,"max_iter":100, "MaxIter": 100, "MaxIt" : 100, "print_level":0,"derivative_test":"first-order", "hessian_approximation": "exact", "UserHM": True, "OptTolAbs":1e-10, "OptTol":1e-10}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput([0,1],"x0")
      solver.setInput([-10,1.2],"lbx")
      solver.setInput([10,2],"ubx")
      solver.setInput([-10],"lbg")
      solver.setInput([10],"ubg")
      solver.evaluate()
      if float(solver.getOutput("x")[0])<0: # JOEL: There appears to be two local minima
        self.assertAlmostEqual(solver.getOutput("f")[0],4.3817250416084308,str(Solver))
        self.assertAlmostEqual(solver.getOutput("x")[0],-1.0910624688699295,6,str(Solver))
        self.assertAlmostEqual(solver.getOutput("x")[1],1.2,5,str(Solver))
        self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,5 if "stabilizedsqp"==Solver else 8,str(Solver)+str(solver_options))
        self.assertAlmostEqual(solver.getOutput("lam_x")[1],-1.9165378046901287,4,str(Solver))
        self.assertAlmostEqual(solver.getOutput("lam_g")[0],0,8,str(Solver))
      else:
        self.assertAlmostEqual(solver.getOutput("f")[0],9.0908263002590e-3,6,str(Solver))
        self.assertAlmostEqual(solver.getOutput("x")[0],1.0952466252248,6,str(Solver))
        self.assertAlmostEqual(solver.getOutput("x")[1],1.2,5,str(Solver))
        self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,5 if "stabilizedsqp"==Solver else 8,str(Solver)+str(solver_options))
        self.assertAlmostEqual(solver.getOutput("lam_x")[1],-8.6963632695079e-2,4,str(Solver))
        self.assertAlmostEqual(solver.getOutput("lam_g")[0],0,8,str(Solver))

  def testactiveLBG(self):
    self.message("active LBG")
    x=SX.sym("x")
    y=SX.sym("y")
    
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=(1-x)**2+100*(y-x**2)**2,g=x+y))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
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
      solver.evaluate()
      self.assertAlmostEqual(solver.getOutput("f")[0],4.252906468284e-3,6,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[0],1.065181061847138,6,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[1],1.1348189166291160,6,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,8,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,4,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_g")[0],-4.1644422845712e-2,3,str(Solver))

  def testactiveUBG(self):
    self.message("active UBG")
    x=SX.sym("x")
    y=SX.sym("y")
    
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=(1-x)**2+100*(y-x**2)**2,g=x+y))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
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
      solver.evaluate()
      self.assertAlmostEqual(solver.getOutput("f")[0],4.64801220074552e-3,6,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[0],9.318651964592811e-1,5,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[1],8.68134821123689e-1,5,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,8,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,4,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_g")[0],4.75846495145007e-2,5,str(Solver))
      
  def testactiveUBX(self):
    self.message("active UBX")
    x=SX.sym("x")
    y=SX.sym("y")
    
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=(1-x)**2+100*(y-x**2)**2,g=x+y))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
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
      solver.evaluate()
      self.assertAlmostEqual(solver.getOutput("f")[0],2.626109721583e-3,6,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[0],9.4882542279172277e-01,6,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[1],0.9,6,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,8,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],5.39346608659e-2,4,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_g")[0],0,8,str(Solver))
      
  def test_QP(self):
    self.message("QP")

    N = 50 

    x = SX.sym("x",N)
    x0 = DMatrix(range(N))
    H = diag(range(1,N+1))
    obj = 0.5*mul([(x-x0).T,H,(x-x0)])

    nlp = SXFunction(nlpIn(x=x),nlpOut(f=obj))
    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-8,"tol_pr":1e-10,"TolOpti":1e-25,"hessian_approximation":"limited-memory","max_iter":100, "MaxIter": 100,"print_level":0}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
      solver.init()
      solver.setInput(-1000,"lbx")
      solver.setInput(1000,"ubx")
      solver.evaluate()
      self.checkarray(solver.getOutput("x"),x0,str(Solver),digits=2)
      self.assertAlmostEqual(solver.getOutput("f")[0],0,3,str(Solver))
      self.checkarray(solver.getOutput("lam_x"),DMatrix.zeros(N,1),str(Solver),digits=4)
      
      
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

    x=SX.sym("x",2)
    nlp=SXFunction(nlpIn(x=x),nlpOut(f=0.5*mul([x.T,H,x])+mul(G.T,x),g=mul(A,x)))

    for Solver, solver_options in solvers:
      self.message(str(Solver))
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-8,"tol_pr":1e-10,"TolOpti":1e-25,"hessian_approximation":"limited-memory","max_iter":100,"MaxIter": 100,"print_level":0, "fixed_variable_treatment": "make_constraint"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
          
      solver.init()
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lbg")
      solver.setInput(UBA,"ubg")

      solver.evaluate()

      self.assertAlmostEqual(solver.getOutput()[0],0.5,6,str(Solver))
      self.assertAlmostEqual(solver.getOutput()[1],1.25,6,str(Solver))
    
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],4.75,6,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,6,str(Solver))

      self.checkarray(solver.getOutput("lam_g"),DMatrix([0,2,0]),str(Solver),digits=6)
      
      self.assertAlmostEqual(solver.getOutput("f")[0],-7.4375,6,str(Solver))
      
  def test_QP2(self):
    H = DMatrix([[1,-1],[-1,2]])
    G = DMatrix([-2,-6])
    A =  DMatrix([[1, 1],[-1, 2],[2, 1]])

    LBA = DMatrix([-inf]*3)
    UBA = DMatrix([2, 2, 3])

    LBX = DMatrix([0.5,0])
    UBX = DMatrix([0.5,inf])

    x=SX.sym("x",2)
    nlp=SXFunction(nlpIn(x=x),nlpOut(f=0.5*mul([x.T,H,x])+mul(G.T,x),g=mul(A,x)))

    for Solver, solver_options in solvers:
      self.message(Solver)
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-8,"TolOpti":1e-25,"hessian_approximation":"limited-memory","max_iter":100,"MaxIter": 100,"print_level":0, "fixed_variable_treatment": "make_constraint"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
          
      solver.init()
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lbg")
      solver.setInput(UBA,"ubg")
      if 'sqic' in str(solver_options):
        continue
      if Solver=='worhp':
        with self.assertRaises(Exception):
          solver.evaluate()
        return

      solver.evaluate()

      self.assertAlmostEqual(solver.getOutput("x")[0],0.5,6,str(Solver))
      self.assertAlmostEqual(solver.getOutput("x")[1],1.25,6,str(Solver))
    
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],4.75,6,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,6,str(Solver))

      self.checkarray(solver.getOutput("lam_g"),DMatrix([0,2,0]),str(Solver),digits=6)
      
      self.assertAlmostEqual(solver.getOutput("f")[0],-7.4375,6,str(Solver))
      
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-8,"TolOpti":1e-25,"hessian_approximation":"exact","UserHM":True,"max_iter":100,"MaxIter": 100,"print_level":0, "fixed_variable_treatment": "make_constraint"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
          
      solver.init()
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lbg")
      solver.setInput(UBA,"ubg")

      solver.evaluate()

      self.assertAlmostEqual(solver.getOutput()[0],0.5,6,str(Solver))
      self.assertAlmostEqual(solver.getOutput()[1],1.25,6,str(Solver))
    
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],4.75,6,str(Solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,6,str(Solver))

      self.checkarray(solver.getOutput("lam_g"),DMatrix([0,2,0]),str(Solver),digits=6)
      
      self.assertAlmostEqual(solver.getOutput("f")[0],-7.4375,6,str(Solver))

  def test_QP2_unconvex(self):
    H = DMatrix([[1,-1],[-1,-2]])
    G = DMatrix([-2,-6])
    A =  DMatrix([[1, 1],[-1, 2],[2, 1]])
    
    LBA = DMatrix([-inf]*3)
    UBA = DMatrix([2, 2, 3])

    LBX = DMatrix([0]*2)
    UBX = DMatrix([inf]*2)

    x=SX.sym("x",2)
    nlp=SXFunction(nlpIn(x=x),nlpOut(f=0.5*mul([x.T,H,x])+mul(G.T,x),g=mul(A,x)))

    for Solver, solver_options in solvers:
      self.message(Solver)
      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-8,"TolOpti":1e-25,"hessian_approximation":"limited-memory","max_iter":100,"MaxIter": 100,"print_level":0, "fixed_variable_treatment": "make_constraint"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
          
      solver.init()
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lbg")
      solver.setInput(UBA,"ubg")

      solver.evaluate()

      self.assertAlmostEqual(solver.getOutput("x")[0],2.0/3,6,str(solver))
      self.assertAlmostEqual(solver.getOutput("x")[1],4.0/3,6,str(solver))
    
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,6,str(solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,6,str(solver))

      self.checkarray(solver.getOutput("lam_g"),DMatrix([4+8.0/9,20.0/9,0]),str(solver),digits=6)
      
      self.assertAlmostEqual(solver.getOutput("f")[0],-10-16.0/9,6,str(solver))

      solver = NlpSolver(Solver, nlp)
      solver.setOption(solver_options)
      for k,v in ({"tol":1e-8,"TolOpti":1e-25,"hessian_approximation":"exact","UserHM":True,"max_iter":100,"MaxIter": 100,"print_level":0, "fixed_variable_treatment": "make_constraint"}).iteritems():
        if solver.hasOption(k):
          solver.setOption(k,v)
          
      solver.init()
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lbg")
      solver.setInput(UBA,"ubg")

      solver.evaluate()

      self.assertAlmostEqual(solver.getOutput()[0],2.0/3,6,str(solver))
      self.assertAlmostEqual(solver.getOutput()[1],4.0/3,6,str(solver))
    
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,6,str(solver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,6,str(solver))

      self.checkarray(solver.getOutput("lam_g"),DMatrix([4+8.0/9,20.0/9,0]),str(solver),digits=6)
      
      self.assertAlmostEqual(solver.getOutput("f")[0],-10-16.0/9,6,str(solver))
      
  def test_bug(self):
    x = MX.sym("x", 3)
    y = MX.sym("y", 2)
    f = MXFunction([x, y], [1.])
    f.init()
    aa = MX.sym("aa", 5)
    a = aa[:3]
    b = aa[3:]
    [f_call] = f.call([a, b])
    nlp = MXFunction(nlpIn(x=aa), nlpOut(f=f_call))
    for Solver, solver_options in solvers:
      solver = NlpSolver(Solver, nlp)
      solver = NlpSolver("ipopt", nlp)
      solver.init() 
      
  @requiresPlugin(NlpSolver,"snopt")
  def test_permute(self):
    for Solver, solver_options in solvers:
      if "snopt" not in str(Solver): continue
      for permute_g in itertools.permutations(range(3)):
        for permute_x in itertools.permutations(range(4)):
          x=SX.sym("x",4)
          x1,x2,x3,x4 = x[permute_x]
          g = [x1**2+x2**2+x3,
              x2**4+x4,
              2*x1+4*x2]
          f= (x1+x2+x3)**2+3*x3+5*x4
          F= SXFunction(nlpIn(x=x),nlpOut(f=f,g=vertcat(g)[permute_g]))
          F.init()
          
          solver = NlpSolver(Solver,F)
          solver.setOption(solver_options)
          solver.init()
          
          
          ubx = solver.getInput("ubx")
          ubx[permute_x]= DMatrix([inf,inf,inf,inf])
          solver.setInput(ubx,"ubx")
          
          
          lbx = solver.getInput("lbx")
          lbx[permute_x]= DMatrix([-inf,-inf,0,0])
          solver.setInput(lbx,"lbx")
          
          solver.setInput(DMatrix([2,4,inf])[permute_g],"ubg")
          solver.setInput(DMatrix([2,4,0])[permute_g],"lbg")
          
          x0 = solver.getInput("x0")
          x0[permute_x] = DMatrix([-0.070,1.41,0,0.0199])
          solver.setInput(x0,"x0")
          
          solver.evaluate()


          F.setInput(solver.getOutput("x"))
          F.evaluate()

          self.checkarray(solver.getOutput("f"),DMatrix([1.9001249992187681e+00]),digits=7)
          self.checkarray(solver.getOutput("x")[permute_x],DMatrix([-7.0622015054877127e-02,1.4124491251068008e+00,0,1.9925001159906402e-02]),failmessage=str(permute_x)+str(permute_g),digits=7)
          self.checkarray(solver.getOutput("lam_x")[permute_x],DMatrix([0,0,-2.4683779218120115e+01,0]),digits=7)
          self.checkarray(solver.getOutput("lam_g"),DMatrix([1.9000124997534527e+01,-5,0])[permute_g],digits=7)
          self.checkarray(solver.getOutput("g"),DMatrix([2,4,5.5085524702939])[permute_g],digits=7)

  @requiresPlugin(NlpSolver,"snopt")
  def test_permute2(self):
    for Solver, solver_options in solvers:
      if "snopt" not in str(Solver): continue
      for permute_g in itertools.permutations(range(3)):
        for permute_x in itertools.permutations(range(4)):
          x=SX.sym("x",4)
          x1,x2,x3,x4 = x[permute_x]
          g = [x1**2+x2+x3,
              x3**2+x4,
              2*x1+4*x2]
          f= x1**2+x3**2
          F= SXFunction(nlpIn(x=x),nlpOut(f=f,g=vertcat(g)[permute_g]))
          F.init()
          
          solver = NlpSolver(Solver,F)
          solver.setOption(solver_options)
          solver.init()

          ubx = solver.getInput("ubx")
          ubx[permute_x]= DMatrix([inf,inf,inf,inf])
          solver.setInput(ubx,"ubx")
          
          lbx = solver.getInput("lbx")
          lbx[permute_x]= DMatrix([-inf,-inf,0,0])
          solver.setInput(lbx,"lbx")

          solver.setInput(DMatrix([2,4,inf])[permute_g],"ubg")
          solver.setInput(DMatrix([2,4,0])[permute_g],"lbg")
          
          x0 = solver.getInput("x0")
          
          x0[permute_x] = DMatrix([-0.070,1.41,0,0.0199])
          solver.setInput(x0,"x0")
          
          solver.evaluate()

          self.checkarray(solver.getOutput("f"),DMatrix([0]),digits=8)
          self.checkarray(solver.getOutput("x")[permute_x],DMatrix([0,2,0,4]),digits=4,failmessage=str(permute_x)+str(permute_g))
          self.checkarray(solver.getOutput("lam_x")[permute_x],DMatrix([0,0,0,0]),digits=3)
          self.checkarray(solver.getOutput("lam_g"),DMatrix([0,0,0])[permute_g],digits=3)
          #self.checkarray(solver.getOutput("g"),DMatrix([2,4,5.50855])[permute_g])

  @requiresPlugin(NlpSolver,"snopt")
  def test_permute3(self):
    for Solver, solver_options in solvers:
      if "snopt" not in str(Solver): continue
      for permute_g in itertools.permutations(range(3)):
        for permute_x in itertools.permutations(range(4)):
          x=SX.sym("x",4)
          x1,x2,x3,x4 = x[permute_x]
          g = [x1**2+x2+x3,
              x3**2+x4,
              2*x1+4*x2]
          f= x1**2+x3**2+2*x2
          F= SXFunction(nlpIn(x=x),nlpOut(f=f,g=vertcat(g)[permute_g]))
          F.init()
          
          solver = NlpSolver(Solver,F)
          solver.setOption(solver_options)
          solver.init()
          
          ubx = solver.getInput("ubx")
          ubx[permute_x]= DMatrix([inf,inf,inf,inf])
          solver.setInput(ubx,"ubx")

          lbx = solver.getInput("lbx")
          lbx[permute_x]= DMatrix([-inf,-inf,0,0])
          solver.setInput(lbx,"lbx")
          
          solver.setInput(DMatrix([2,4,inf])[permute_g],"ubg")
          solver.setInput(DMatrix([2,4,0])[permute_g],"lbg")
          
          x0 = solver.getInput("x0") 
          x0[permute_x] = DMatrix([1,-0.5,0.5,4])
          solver.setInput(x0,"x0")
          
          solver.evaluate()

          self.checkarray(solver.getOutput("f"),DMatrix([9.9030108869944522e-01]),failmessage=str(permute_x)+str(permute_g))
          self.checkarray(solver.getOutput("x")[permute_x],DMatrix([1.53822842722,-0.76911421361,0.402967519303,3.83761717839]),digits=6)
          self.checkarray(solver.getOutput("lam_x")[permute_x],DMatrix([0,0,0,0]),digits=7)
          self.checkarray(solver.getOutput("lam_g"),DMatrix([-8.0593503860219973e-01,6.52750754744e-10,-0.298516240384])[permute_g],failmessage=str(permute_x)+str(permute_g),digits=8)
          #self.checkarray(solver.getOutput("g"),DMatrix([2,4,5.50855])[permute_g])
        
  @requiresPlugin(NlpSolver,"snopt")
  def test_classifications(self):      
    x=SX.sym("x")
    y=SX.sym("y")
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=(1-x)**2+7.7*y,g=y**2))

    solver = NlpSolver("snopt", nlp)
        
    #solver.setOption("detect_linear",False)
    solver.setOption("verbose",True)
    solver.setOption("monitor",["setup_nlp","eval_nlp"])
    solver.setOption("Verify level",3)
    #solver.setOption("_iteration_limit",1000)

    solver.init()
    solver.setInput([1,1],"x0")
    solver.setInput([-10,0],"lbx")
    solver.setInput([10,2],"ubx")
    solver.setInput([-10],"lbg")
    solver.setInput([10],"ubg")

    solver.evaluate()
    
    self.checkarray(solver.getOutput("f"),DMatrix([0]))
    self.checkarray(solver.getOutput("x"),DMatrix([1,0]))
    self.checkarray(solver.getOutput("lam_x"),DMatrix([0,-7.7]),digits=7)
    self.checkarray(solver.getOutput("lam_g"),DMatrix([0]))
    
  def test_pathological(self):      
    x=SX.sym("x")
    y=SX.sym("y")
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=(1-x)**2+y**2))

    for Solver, solver_options in solvers:
      self.message(str(Solver))
      if "worhp"==Solver or "stabilizedsqp"==Solver : continue
      solver = NlpSolver(Solver, nlp)
          
      #solver.setOption("detect_linear",False)
      solver.setOption("verbose",True)
      solver.setOption(solver_options)
      
      solver.init()
      solver.setInput([1,1],"x0")
      solver.setInput([-10,-1],"lbx")
      solver.setInput([10,2],"ubx")

      solver.evaluate()
      
      self.checkarray(solver.getOutput("f"),DMatrix([0]),digits=7)
      self.checkarray(solver.getOutput("x"),DMatrix([1,0]),digits=7,failmessage=str(Solver))
      self.checkarray(solver.getOutput("lam_x"),DMatrix([0,-0]),digits=7,failmessage=str(Solver))

  def test_pathological2(self):      
    x=SX.sym("x")
    y=SX.sym("y")
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=(1-x)**2+y))

    for Solver, solver_options in solvers:
      self.message(Solver)
      solver = NlpSolver(Solver, nlp)
          
      #solver.setOption("detect_linear",False)
      solver.setOption("verbose",True)
      solver.setOption(solver_options)

      solver.init()
      solver.setInput([1,1],"x0")
      solver.setInput([-10,0],"lbx")
      solver.setInput([10,2],"ubx")

      solver.evaluate()
      
      self.checkarray(solver.getOutput("f"),DMatrix([0]),digits=7)
      self.checkarray(solver.getOutput("x"),DMatrix([1,0]),digits=7)
      self.checkarray(solver.getOutput("lam_x"),DMatrix([0,-1]),digits=7)

  def test_pathological3(self):      
    x=SX.sym("x")
    y=SX.sym("y")
    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=(1-x)**2,g=x+y))

    for Solver, solver_options in solvers:
      self.message(str(Solver))
      if "worhp"==Solver: continue
      solver = NlpSolver(Solver, nlp)
          
      #solver.setOption("detect_linear",False)
      solver.setOption("verbose",True)
      solver.setOption(solver_options)

      solver.init()
      solver.setInput([1,1],"x0")
      solver.setInput([-10,0],"lbx")
      solver.setInput([10,2],"ubx")
      solver.setInput([2],"lbg")
      solver.setInput([2],"ubg")
      
      solver.evaluate()
      
      self.checkarray(solver.getOutput("f"),DMatrix([0]),digits=7)
      self.checkarray(solver.getOutput("x"),DMatrix([1,1]),digits=7)
      self.checkarray(solver.getOutput("lam_x"),DMatrix([0,0]),digits=7)
    
  def test_pathological4(self):      
    x=SX.sym("x")
    nlp=SXFunction(nlpIn(x=x),nlpOut(f=x*x))

    for Solver, solver_options in solvers:
      self.message(Solver)
      if "worhp"==Solver: continue
      solver = NlpSolver(Solver, nlp)
          
      #solver.setOption("detect_linear",False)
      solver.setOption("verbose",True)
      solver.setOption(solver_options)

      solver.init()
      solver.setInput([0],"x0")
      solver.setInput([0],"lbx")
      solver.setInput([0],"ubx")
      
      solver.evaluate()
      
      self.checkarray(solver.getOutput("f"),DMatrix([0]),digits=7)
      self.checkarray(solver.getOutput("x"),DMatrix([0]),digits=7)
      self.checkarray(solver.getOutput("lam_x"),DMatrix([0]),digits=7)
      
if __name__ == '__main__':
    unittest.main()
    print solvers

