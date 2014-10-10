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
import numpy as n
import unittest
from types import *
from helpers import *
import copy

scipy_available = True
try:
	import scipy.special
	from scipy.linalg import expm
except:
	scipy_available = False
	
integrators = []

try:
  Integrator.loadPlugin("cvodes")
  integrators.append(("cvodes",["ode"],{"abstol": 1e-15,"reltol":1e-15,"fsens_err_con": True,"quad_err_con": False}))
except:
  pass
  
try:
  Integrator.loadPlugin("idas")
  integrators.append(("idas",["dae","ode"],{"abstol": 1e-15,"reltol":1e-15,"fsens_err_con": True,"calc_icB":True}))
except:
  pass

integrators.append(("collocation",["dae","ode"],{"implicit_solver":"kinsol","number_of_finite_elements": 18}))

try:
  Integrator.loadPlugin("oldcollocation")
  integrators.append(("oldcollocation",["dae","ode"],{"implicit_solver":"kinsol","number_of_finite_elements": 18,"startup_integrator":"cvodes"}))
  #integrators.append(("oldcollocation",["dae","ode"],{"implicit_solver":"nlp","number_of_finite_elements": 100,"startup_integrator":"cvodes","implicit_solver_options": {"nlp_solver": "ipopt","linear_solver_creator": "csparse"}}))
except:
  pass

integrators.append(("rk",["ode"],{"number_of_finite_elements": 1000}))

print "Will test these integrators:"
for cl, t, options in integrators:
  print cl, " : ", t

class Integrationtests(casadiTestCase):

  @slow()
  def test_full(self):
    num = self.num
    tc = DMatrix(n.linspace(0,num['tend'],100))
    
    t=SX.sym("t")
    q=SX.sym("q")
    p=SX.sym("p")
    
    out = SXFunction(daeIn(t=t, x=q, p=p),[q,t,p])
    out.init()
        
    f=SXFunction(daeIn(t=t, x=q, p=p),daeOut(ode=q/p*t**2))
    f.init()
    integrator = Integrator("cvodes",f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    integrator.setOption("fsens_err_con", True)
    #integrator.setOption("verbose",True)
    integrator.setOption("t0",0)
    integrator.setOption("tf",2.3)
    integrator.init()
    tf = 2.3
    
    solution = SXFunction(integratorIn(x0=q, p=p),integratorOut(xf=q*exp(tf**3/(3*p))))
    solution.init()
    
    for f in [solution,integrator]:
      f.setInput(0.3,"x0")
      f.setInput(0.7,"p")
    
    self.checkfunction(integrator,solution,digits=6)

  def test_tools_trivial(self):
    num = self.num

    t=SX.sym("t")
    q=SX.sym("q")
    p=SX.sym("p")
    
    f=SXFunction(daeIn(x=q),daeOut(ode=q))
    f.init()
    tf = 1
    
    for integrator in [
         explicitRK(f,tf,4,10),
         implicitRK(f,"newton",{"linear_solver": "csparse"},tf,4,"radau",10)
       ]:
      integrator.init()
      
      solution = SXFunction(integratorIn(x0=q),integratorOut(xf=q*exp(tf)))
      solution.init()
      
      for f in [solution,integrator]:
        f.setInput(1,"x0")
        
      integrator.evaluate()
      

      self.checkfunction(integrator,solution,digits=5)

  @slow()
  def test_tools(self):
    num = self.num

    t=SX.sym("t")
    q=SX.sym("q")
    p=SX.sym("p")
    
    out = SXFunction(daeIn(t=t, x=q, p=p),[q,t,p])
    out.init()

    f=SXFunction(daeIn(t=t, x=q, p=p),daeOut(ode=q/p*t**2))
    f.init()
    tf = 1
    for integrator in [
         explicitRK(f,tf,4,500),
         implicitRK(f,"newton",{"linear_solver": "csparse"},tf,4,"radau",50)
       ]:
      integrator.init()
      
      solution = SXFunction(integratorIn(x0=q, p=p),integratorOut(xf=q*exp(tf**3/(3*p))))
      solution.init()
      
      for f in [solution,integrator]:
        f.setInput(0.3,"x0")
        f.setInput(0.7,"p")
      
      self.checkfunction(integrator,solution,digits=4)
    
  @memory_heavy()
  def test_jac(self):
    self.message("Test exact jacobian #536")
    # This test is not automized, but works by inspection only.
    # To activate, recompile after ucnommenting the printout lines in cvodes.c, near "Used for validating casadi#536"
    #return
    DMatrix.setPrecision(18)

    tstart = SX.sym("tstart")
    tend = SX.sym("tend")
    
    integrators = [
              ("idas",["dae","ode"],{"abstol": 1e-9,"reltol":1e-9,"fsens_err_con": True,"calc_ic":True,"calc_icB":True}),
              ("cvodes",["ode"],{"abstol": 1e-5,"reltol":1e-5,"fsens_err_con": False,"quad_err_con": False})
              ]

    def variations(p_features, din, dout, rdin, rdout, *args):
      if "ode" in p_features:
        p_features_ = copy.copy(p_features)
        p_features_[p_features.index("ode")] = "dae"
        din_ = copy.copy(din)
        dout_ = copy.copy(dout)
        rdin_ = copy.copy(rdin)
        rdout_ = copy.copy(rdout)
        z = SX.sym("x", din_["x"].shape)
        din_["z"] = z
        dout_["ode"] = z
        dout_["alg"] = ( dout["ode"] - z) * (-0.8)
        if len(rdin_)>0:
          rz = SX.sym("rx", rdin_["rx"].shape)
          rdin_["rz"] = rz
          rdin_["z"] = z
          rdout_["ode"] = rz
          rdout_["alg"] = ( rdout["ode"] - rz) * (-0.7)
          
        yield (p_features, din, dout, rdin, rdout) + tuple(args)
        yield (p_features_, din_, dout_, rdin_, rdout_) + tuple(args)
      else:
        yield (p_features, din, dout, rdin, rdout) + tuple(args)
        
    def checks(): 
      Ns = 1
      
      x  = SX.sym("x")
      rx = SX.sym("rx")
      t = SX.sym("t")

      ti = (0,0.9995)
      pointA = {'x0': 1, 'rx0': 1}
      
      si = {'x0':x, 'rx0': rx}
      
      #sol = {'rxf': 1.0/(1-tend)}
      sol = {'rxf': rx*exp(tend), 'xf': x*exp(tend)}
     
      yield (["ode"],{'x':x,'t':t},{'ode':x},{'x':x,'rx':rx,'t':t},{'ode': rx},si,sol,pointA,ti)
      
      
    refXF = refRXF = None

    for tt in checks():
      for p_features, din, dout, rdin, rdout,  solutionin, solution, point, (tstart_, tend_) in variations(*tt):
        for Integrator, features, options in integrators:
          self.message(Integrator)
          dummyIntegrator = c.Integrator(Integrator,c.SXFunction())
          if p_features[0] in features:
            g = Function()
            if len(rdin)>1:
              g = SXFunction(rdaeIn(**rdin),rdaeOut(**rdout))
              g.init()
               
            f = SXFunction(daeIn(**din),daeOut(**dout))
            f.init()
            
            for k in solution.keys():
              solution[k] = substitute(solution[k],vertcat([tstart,tend]),vertcat([tstart_,tend_]))

            fs = SXFunction(integratorIn(**solutionin),integratorOut(**solution))
            fs.init()
              
          
            def itoptions(post=""):
              yield {"iterative_solver"+post: "gmres"}
              yield {"iterative_solver"+post: "bcgstab"}
              yield {"iterative_solver"+post: "tfqmr", "use_preconditionerB": True, "linear_solverB" : "csparse"} # Bug in Sundials? Preconditioning seems to be needed
             
            def solveroptions(post=""):
              yield {"linear_solver_type" +post: "dense" }
              allowedOpts = list(dummyIntegrator.getOptionAllowed("linear_solver_type" +post))
              #allowedOpts.remove("iterative") # disabled, see #1231
              if "iterative" in allowedOpts:
                  for it in itoptions(post):
                      d = {"linear_solver_type" +post: "iterative" }
                      d.update(it)
                      yield d
              if "banded" in allowedOpts:
                  yield {"linear_solver_type" +post: "banded" }
              yield {"linear_solver_type" +post: "user_defined", "linear_solver"+post: "csparse" }
                
            for a_options in solveroptions("B"):
              for f_options in solveroptions():
                message = "f_options: %s , a_options: %s" % (str(f_options) , str(a_options))
                print message
                integrator = c.Integrator(Integrator,f,g)
                integrator.setOption("exact_jacobianB",True)
                integrator.setOption("gather_stats",True)
                #integrator.setOption("verbose",True)
                #integrator.setOption("monitor",["djacB","resB","djac","res"])
                integrator.setOption("t0",tstart_)
                integrator.setOption("tf",tend_)
                integrator.setOption(options)
                integrator.setOption(f_options)
                integrator.setOption(a_options)
                integrator.init()
                for ff in [fs,integrator]:
                  for k,v in point.items():
                    i = getattr(casadi,('integrator_'+k).upper())
                    if not ff.getInput(i).isEmpty():
                      ff.setInput(v,i)

                integrator.evaluate()
                fs.evaluate()
                print "res=",integrator.getOutput("xf")-fs.getOutput("xf"), fs.getOutput("xf")
                print "Rres=",integrator.getOutput("rxf")-fs.getOutput("rxf"), fs.getOutput("rxf")
                # self.checkarray(integrator.getOutput("rxf"),fs.getOutput("rxf"),digits=4)
                stats = integrator.getStats()
                
                print stats
                self.assertTrue(stats["nsteps"]<1500)
                self.assertTrue(stats["nstepsB"]<2500)
                self.assertTrue(stats["nlinsetups"]<100)
                self.assertTrue(stats["nlinsetupsB"]<250)

  @memory_heavy()
  def test_lsolvers(self):
    self.message("Test different linear solvers")

    tstart = SX.sym("tstart")
    tend = SX.sym("tend")
    
    integrators = [
              ("idas",["dae","ode"],{"abstol": 1e-9,"reltol":1e-9,"fsens_err_con": True,"calc_ic":True,"calc_icB":True}),
              ("cvodes",["ode"],{"abstol": 1e-15,"reltol":1e-15,"fsens_err_con": True,"quad_err_con": False})
              ]
              
    def checks():  
      t=SXElement.sym("t")
      x=SXElement.sym("x")
      rx=SXElement.sym("rx")
      p=SXElement.sym("p")
      dp=SXElement.sym("dp")

      z=SXElement.sym("z")
      rz=SXElement.sym("rz")
      rp=SXElement.sym("rp")    
      solutionin = {'x0':x, 'p': p, 'rx0': rx,'rp' : rp}            
      pointA = {'x0':7.1,'p': 2, 'rx0': 0.13, 'rp': 0.127}
      ti = (0.2,2.3)
      yield (["dae"],{'x': x, 'z': z},{'alg': x-z, 'ode': z},{'x': x, 'z': z, 'rx': rx, 'rz': rz},{'alg': x-rz, 'ode': rz},solutionin,{'rxf': rx+x*(exp(tend-tstart)-1), 'xf':x*exp(tend-tstart)},pointA,ti)
      if not(args.run_slow): return
      yield (["dae"],{'x': x, 'z': z},{'alg': x-z, 'ode': z},{'x': x, 'z': z, 'rx': rx, 'rz': rz},{'alg': rx-rz, 'ode': rz},solutionin,{'rxf': rx*exp(tend-tstart), 'xf':x*exp(tend-tstart)},pointA,ti)
      yield (["ode"],{'x': x},{'ode': x},{'x': x,'rx': rx},{'ode': x},solutionin,{'rxf': rx+x*(exp(tend-tstart)-1), 'xf':x*exp(tend-tstart)},pointA,ti)
      yield (["ode"],{'x': x},{'ode': x},{'x': x,'rx': rx},{'ode': rx},solutionin,{'rxf': rx*exp(tend-tstart), 'xf':x*exp(tend-tstart)},pointA,ti)
      
      A=array([1,0.1])
      p0 = 1.13

      q=SX.sym("y",2,1)
      y0=q[0]
      yc0=dy0=q[1]
      p=SX.sym("p",1,1)

      s1=(2*y0-log(yc0**2/p+1))/2-log(cos(arctan(yc0/sqrt(p))+sqrt(p)*(tend-tstart)))
      s2=sqrt(p)*tan(arctan(yc0/sqrt(p))+sqrt(p)*(tend-tstart))
      yield (["ode"],{'x':q,'p':p},{'ode': vertcat([q[1],p[0]+q[1]**2 ])},{},{},{'x0':q, 'p': p} ,{'xf': vertcat([s1,s2])},{'x0': A, 'p': p0},(0,0.4) )

    for p_features, din, dout, rdin, rdout, solutionin, solution, point, (tstart_, tend_) in checks():

      for Integrator, features, options in integrators:
        self.message(Integrator)
        dummyIntegrator = c.Integrator(Integrator,SXFunction())
        if p_features[0] in features:
          g = Function()
          if len(rdin)>1:
            g = SXFunction(rdaeIn(**rdin),rdaeOut(**rdout))
            g.init()
             
          f = SXFunction(daeIn(**din),daeOut(**dout))
          f.init()
            
          for k in solution.keys():
            solution[k] = substitute(solution[k],vertcat([tstart,tend]),vertcat([tstart_,tend_]))
          
          fs = SXFunction(integratorIn(**solutionin),integratorOut(**solution))
          fs.init()
        
          def itoptions(post=""):
            yield {"iterative_solver"+post: "gmres"}
            yield {"iterative_solver"+post: "bcgstab"}
            yield {"iterative_solver"+post: "tfqmr", "use_preconditionerB": True, "linear_solverB" : "csparse"} # Bug in Sundials? Preconditioning seems to be needed
           
          def solveroptions(post=""):
            yield {"linear_solver_type" +post: "dense" }
            allowedOpts = list(dummyIntegrator.getOptionAllowed("linear_solver_type" +post))
            #allowedOpts.remove("iterative")  # disabled, see #1231
            if "iterative" in allowedOpts:
                for it in itoptions(post):
                    d = {"linear_solver_type" +post: "iterative" }
                    d.update(it)
                    yield d
            if "banded" in allowedOpts:
                yield {"linear_solver_type" +post: "banded" }
            yield {"linear_solver_type" +post: "user_defined", "linear_solver"+post: "csparse" }
              
          for a_options in solveroptions("B"):
            for f_options in solveroptions():
              message = "f_options: %s , a_options: %s" % (str(f_options) , str(a_options))
              print message
              integrator = c.Integrator(Integrator,f,g)
              integrator.setOption("exact_jacobianB",True)
              integrator.setOption("t0",tstart_)
              integrator.setOption("tf",tend_)
              integrator.setOption(options)
              integrator.setOption(f_options)
              integrator.setOption(a_options)
              integrator.init()
              
              for ff in [fs,integrator]:
                for k,v in point.items():
                  i = getattr(casadi,('integrator_'+k).upper())
                  if not ff.input(i).isEmpty():
                    ff.setInput(v,i)

              integrator.evaluate()
              
              self.checkfunction(integrator,fs,gradient=False,hessian=False,sens_der=False,evals=False,digits=4,digits_sens=4,failmessage=message,verbose=False)
              
              


  @memory_heavy()
  def test_X(self):
    self.message("Extensive integrator tests")
    
    num=self.num
    tstart = SX.sym("tstart")
    tend = SX.sym("tstart")

    
    for Integrator, features, options in integrators:
      self.message(Integrator)
        
        
      def variations(p_features, din, dout, rdin, rdout, *args):
        if "ode" in p_features:
          p_features_ = copy.copy(p_features)
          p_features_[p_features.index("ode")] = "dae"
          din_ = copy.copy(din)
          dout_ = copy.copy(dout)
          rdin_ = copy.copy(rdin)
          rdout_ = copy.copy(rdout)
          z = SX.sym("x", din_["x"].shape)
          din_["z"] = z
          dout_["ode"] = z
          dout_["alg"] = ( dout["ode"] - z) * (-0.8)
          if len(rdin_)>0:
            rz = SX.sym("rx", rdin_["rx"].shape)
            rdin_["rz"] = rz
            rdin_["z"] = z
            rdout_["ode"] = rz
            rdout_["alg"] = ( rdout["ode"] - rz) * (-0.7)
            
          yield (p_features, din, dout, rdin, rdout) + tuple(args)
          yield (p_features_, din_, dout_, rdin_, rdout_) + tuple(args)
        else:
          yield (p_features, din, dout, rdin, rdout) + tuple(args)
        
      def checks():
        x0=num['q0']
        p_=num['p']
        rx0_= 0.13
        t=SX.sym("t")
        x=SX.sym("x")
        rx=SX.sym("rx")
        p=SX.sym("p")
        dp=SX.sym("dp")

        z=SX.sym("z")
        rz=SX.sym("rz")
        rp=SX.sym("rp")
        
        si = {'x0':x, 'p': p, 'rx0': rx,'rp' : rp}            
        pointA = {'x0':x0,'p': p_, 'rx0': rx0_, 'rp': 0.127}
        
        ti = (0.2,num['tend'])
        yield (["ode"],{'x':x},{'ode': 0},{},{},si,{'xf':x},pointA,ti)
        yield (["ode"],{'x':x},{'ode': 1},{},{},si,{'xf':x+(tend-tstart)},pointA,ti)
        yield (["ode"],{'x':x},{'ode': x},{},{},si,{'xf':x*exp(tend-tstart)},pointA,ti)
        yield (["ode"],{'x':x,'t':t},{'ode': t},{},{},si,{'xf':x+(tend**2/2-tstart**2/2)},pointA,ti)
        yield (["ode"],{'x':x,'t':t},{'ode': x*t},{},{},si,{'xf':x*exp(tend**2/2-tstart**2/2)},pointA,ti)
        yield (["ode"],{'x':x,'p':p},{'ode': x/p},{},{},si,{'xf':x*exp((tend-tstart)/p)},pointA,ti)
        if not(args.run_slow): return
        yield (["ode"],{'x':x},{'ode': x,'quad':0},{},{},si,{'qf':0},pointA,ti)
        yield (["ode"],{'x':x},{'ode': x,'quad':1},{},{},si,{'qf':(tend-tstart)},pointA,ti)
        yield (["ode"],{'x':x},{'ode': 0,'quad':x},{},{},si,{'qf':x*(tend-tstart)},pointA,ti)
        #yield ({'x':x},{'ode': 1,'quad':x},{'qf':(x-tstart)*(tend-tstart)+(tend**2/2-tstart**2/2)}), # bug in cvodes quad_err_con
        yield (["ode"],{'x':x},{'ode': x,'quad':x},{},{},si,{'qf':x*(exp(tend-tstart)-1)},pointA,ti)
        yield (["ode"],{'x':x,'t':t},{'ode': x,'quad':t},{},{},si,{'qf':(tend**2/2-tstart**2/2)},pointA,ti)
        yield (["ode"],{'x':x,'t':t},{'ode': x,'quad':x*t},{},{},si,{'qf':x*(exp(tend-tstart)*(tend-1)-(tstart-1))},pointA,ti)
        yield (["ode"],{'x':x,'p':p},{'ode': x,'quad':x/p},{},{},si,{'qf':x*(exp((tend-tstart))-1)/p},pointA,ti)
        yield (["ode"],{'x':x},{'ode':x},{'x':x,'rx':rx},{'ode':0},si,{'rxf': rx},pointA,ti)
        yield (["ode"],{'x':x},{'ode':x},{'x':x,'rx':rx},{'ode':1},si,{'rxf': rx+tend-tstart},pointA,ti)
        yield (["ode"],{'x':x,'t':t},{'ode':x},{'x':x,'rx':rx,'t':t},{'ode':t},si,{'rxf': rx+tend**2/2-tstart**2/2},pointA,ti)
        yield (["ode"],{'x':x},{'ode':x},{'x':x,'rx':rx},{'ode':rx},si,{'rxf': rx*exp(tend-tstart)},pointA,ti)
        yield (["ode"],{'x':x},{'ode':x},{'x':x,'rx':rx},{'ode':x},si,{'rxf': rx+x*(exp(tend-tstart)-1)},pointA,ti)
        yield (["ode"],{'x':x,'t':t},{'ode':x},{'x':x,'rx':rx,'t':t},{'ode':x*t},si,{'rxf': rx+x*(exp(tend-tstart)*(tend-1)-(tstart-1))},pointA,ti)
        yield (["ode"],{'x':x,'t':t},{'ode':x},{'x':x,'rx':rx,'t':t},{'ode':rx*t},si,{'rxf': rx*exp(tend**2/2-tstart**2/2)},pointA,ti)
        yield (["ode"],{'x':x},{'ode':x},{'x':x,'rx':rx},{'ode':rx, 'quad': 0},si,{'rqf': 0},pointA,ti)
        yield (["ode"],{'x':x},{'ode':x},{'x':x,'rx':rx},{'ode':rx, 'quad': 1},si,{'rqf': (tend-tstart)},pointA,ti)
        yield (["ode"],{'x':x},{'ode':x},{'x':x,'rx':rx},{'ode':rx, 'quad': rx},si,{'rqf': rx*(exp(tend-tstart)-1)},pointA,ti)
        yield (["ode"],{'x':x},{'ode':x},{'x':x,'rx':rx},{'ode':rx, 'quad': x},si,{'rqf': x*(exp(tend-tstart)-1)},pointA,ti)
        yield (["ode"],{'x':x,'t':t},{'ode':x},{'x':x,'rx':rx,'t':t},{'ode':rx, 'quad': t},si,{'rqf': (tend**2/2-tstart**2/2)},pointA,ti)
        yield (["ode"],{'x':x,'t':t},{'ode':x},{'x':x,'rx':rx,'t':t},{'ode':rx, 'quad': x*t},si,{'rqf': x*(exp(tend-tstart)*(tend-1)-(tstart-1))},pointA,ti)
        yield (["ode"],{'x':x,'t':t},{'ode':x},{'x':x,'rx':rx,'t':t},{'ode':rx, 'quad': rx*t},si,{'rqf': rx*(exp(tend-tstart)*(tstart+1)-(tend+1))},pointA,ti) # this one is special: integrate(t*rx*exp(tf-t),t,t0,tf)
        yield (["ode"],{'x':x,'p':p},{'ode':x},{'x':x,'rx':rx,'p':p},{'ode':rx, 'quad': p},si,{'rqf': p*(tend-tstart)},pointA,ti)
        yield (["ode"],{'x':x,'p':p},{'ode':x},{'x':x,'rx':rx,'p':p,'rp':rp},{'ode':rx, 'quad': rp},si,{'rqf': rp*(tend-tstart)},pointA,ti)
        yield (["ode"],{'x':x,'t':t},{'ode':x},{'x':x,'rx':rx,'t':t},{'ode':rx*t},si,{'rxf': rx*exp(tend**2/2-tstart**2/2)},pointA,ti)
        yield (["ode"],{'x':x,'t':t},{'ode':x},{'x':x,'rx':rx,'t':t},{'ode':x*t},si,{'rxf': rx+x*(exp(tend-tstart)*(tend-1)-(tstart-1))},pointA,ti)
        yield (["dae"],{'x':x,'z':z},{'ode':z,'alg': -0.8*(z-x),'quad': z},{},{},si,{'qf':x*(exp(tend-tstart)-1)},pointA,ti)
        yield (["dae"],{'x':x,'z':z},{'ode':z,'alg': -0.8*(z-x)},{'x':x,'rx':rx,'rz': rz,'z':z},{'ode':rz, 'alg': -0.7*(rz-rx), 'quad': rz},si,{'rqf': rx*(exp(tend-tstart)-1)},pointA,ti)
        yield (["dae"],{'x':x,'z':z},{'ode':z,'alg': -0.8*(z-x)},{'x':x,'rx':rx,'rz': rz,'z':z},{'ode':rz, 'alg': -0.7*(rz-rx), 'quad': z},si,{'rqf': x*(exp(tend-tstart)-1)},pointA,ti)
        
        
        A=array([1,0.1])
        p0 = 1.13

        q=SX.sym("y",2,1)
        y0=q[0]
        yc0=dy0=q[1]
        p=SX.sym("p",1,1)
        
        s1=(2*y0-log(yc0**2/p+1))/2-log(cos(arctan(yc0/sqrt(p))+sqrt(p)*(tend-tstart)))
        s2=sqrt(p)*tan(arctan(yc0/sqrt(p))+sqrt(p)*(tend-tstart))
        yield (["ode"],{'x':q,'p':p},{'ode': vertcat([q[1],p[0]+q[1]**2 ])},{},{},{'x0':q, 'p': p} ,{'xf': vertcat([s1,s2])},{'x0': A, 'p': p0},(0,0.4) )
      
      for tt in checks():
        print tt
        for p_features, din, dout, rdin, rdout, solutionin, solution, point, (tstart_, tend_) in variations(*tt):
          if p_features[0] in features:
            message = "%s: %s => %s, %s => %s, explicit (%s) tstart = %f" % (Integrator,str(din),str(dout),str(rdin),str(rdout),str(solution),tstart_)
            print message
            g = Function()
            if len(rdin)>1:
              g = SXFunction(rdaeIn(**rdin),rdaeOut(**rdout))
              g.init()
               
            f = SXFunction(daeIn(**din),daeOut(**dout))
            f.init()
            
            for k in solution.keys():
              solution[k] = substitute(solution[k],vertcat([tstart,tend]),vertcat([tstart_,tend_]))
            
            fs = SXFunction(integratorIn(**solutionin),integratorOut(**solution))
            fs.init()
            
            integrator = c.Integrator(Integrator,f,g)
            integrator.setOption(options)
            integrator.setOption("t0",tstart_)
            if integrator.hasOption("abstol"):
              integrator.setOption("abstol",1e-9)
            if integrator.hasOption("reltol"):
              integrator.setOption("reltol",1e-9)
            integrator.setOption("tf",tend_)
            if integrator.hasOption("init_xdot"):
              integrator.setOption("init_xdot",list(DMatrix(point["x0"])))
              integrator.setOption("calc_icB",True)
              integrator.setOption("augmented_options", {"init_xdot":None, "abstol":1e-9,"reltol":1e-9})
            #if "dae" in p_features and integrator.hasOption("init_z"):
            #  integrator.setOption("init_z",[0.1])
            #  integrator.setOption("augmented_options", {"init_z":GenericType(),"init_xdot":GenericType()})
            integrator.init()

#              reproduce = """
#from casadi import *
#t=SXElement.sym("t")
#x=SXElement.sym("x")
#rx=SXElement.sym("rx")
#p=SXElement.sym("p")
#dp=SXElement.sym("dp")

#z=SXElement.sym("z")
#rz=SXElement.sym("rz")
#rp=SXElement.sym("rp")
#f = SXFunction(daeIn(**{din}),daeOut(**{dout}))
#f.init()
#g = SXFunction(rdaeIn(**{rdin}),rdaeOut(**{rdout}))
#g.init()

#integrator = {intclass.__name__}(f,g)
#integrator.setOption({options})
#integrator.init()

#integrator.setInput({x0},"x0")
#if not integrator.input("p").isEmpty():
#  integrator.setInput({p_},"p")
#if not integrator.input("rx0").isEmpty():
#  integrator.setInput(0.13,"rx0")
#if not integrator.input("rp").isEmpty():
#  integrator.setInput(0.127,"rp")
#              """.format(din=din,dout=dout,rdin=rdin,rdout=rdout,x0=x0,p_=p_,intclass=Integrator,options=integrator.dictionary())
#              message+="\nTo reproduce:\n" + reproduce

                
       
            for ff in [fs,integrator]:
              for k,v in point.items():
                i = getattr(casadi,('integrator_'+k).upper())
                if not ff.input(i).isEmpty():
                  ff.setInput(v,i)
            integrator.evaluate()
            
            self.checkfunction(integrator,fs,gradient=False,hessian=False,sens_der=False,evals=False,digits=4,digits_sens=4,failmessage=message,verbose=False)

        
  def setUp(self):
    # Reference solution is x0 e^((t^3-t0^3)/(3 p))
    t=SX.sym("t")
    x=SX.sym("x")
    p=SX.sym("p")
    f=SXFunction(daeIn(t=t, x=x, p=p),daeOut(ode=x/p*t**2))
    f.init()
    integrator = Integrator("cvodes",f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    integrator.setOption("fsens_err_con", True)
    #integrator.setOption("verbose",True)
    integrator.setOption("t0",0)
    integrator.setOption("tf",2.3)
    integrator.init()
    q0   = MX.sym("q0")
    par  = MX.sym("p")
    
    # qend,*_ = integrator.call([q0,par]) # Valid Python3 syntax
    qend, = integratorOut(integrator.call(integratorIn(x0=q0,p=par)),"xf")
    
    qe=MXFunction([q0,par],[qend])
    qe.init()
    self.integrator = integrator
    self.qe=qe
    self.qend=qend
    self.q0=q0
    self.par=par
    self.f = f
    self.num={'tend':2.3,'q0':7.1,'p':2}
    pass
            
  def test_eval2(self):
    self.message('CVodes integration: evaluation with MXFunction indirection')
    num=self.num
    qend=self.qend
    
    par=self.par
    q0=self.q0
    qe=MXFunction([q0,par],[qend[0]])
    qe.init()
    
    f = MXFunction([q0],qe.call([q0,MX(num['p'])]))
    f.init()
    f.setInput([num['q0']],0)
    f.evaluate()

    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(f.getOutput()[0],q0*exp(tend**3/(3*p)),9,"Evaluation output mismatch")
  
  def test_issue92c(self):
    self.message("regression check for issue 92")
    t=SXElement.sym("t")
    x=SXElement.sym("x")
    y=SXElement.sym("y")
    z=x*exp(t)
    f=SXFunction(daeIn(t=t, x=vertcat([x,y])),[vertcat([z,z])])
    f.init()
    # Pass inputs
    f.setInput(1.0,"t")
    f.setInput([1.0,0.0],"x")
    # Evaluate 
    f.evaluate()
    # print result
    print f.getOutput()
  
  def test_issue92b(self):
    self.message("regression check for issue 92")
    t=SXElement.sym("t")
    x=SXElement.sym("x")
    y=SXElement.sym("y")
    f=SXFunction(daeIn(t=t, x=vertcat([x,y])),daeOut(ode=vertcat([x,(1+1e-9)*x])))
    integrator = Integrator("cvodes",f)
    integrator.setOption("fsens_err_con", True)
    integrator.setOption("t0",0)
    integrator.setOption("tf",1)
    integrator.init()
    # Pass inputs
    integrator.setInput([1,0],"x0")
    ## Integrate
    integrator.evaluate()
    # print result
    print integrator.getOutput("xf")
    
  def test_issue92(self):
    self.message("regression check for issue 92")
    t=SXElement.sym("t")
    x=SXElement.sym("x")
    var = MX.sym("var",2,1)

    q = vertcat([x,SXElement.sym("problem")])

    dq=vertcat([x,x])
    f=SXFunction(daeIn(t=t,x=q),daeOut(ode=dq))
    f.init()

    integrator = Integrator("cvodes",f)
    integrator.setOption("fsens_err_con", True)
    integrator.setOption("reltol",1e-12)
    integrator.setOption("t0",0)
    integrator.setOption("tf",1)
    integrator.init()

    qend, = integratorOut(integrator.call(integratorIn(x0=var)),"xf")

    f = MXFunction([var],[qend[0]])
    f.init()

    J=f.jacobian(0)
    J.init()
    J.setInput([1,0])
    J.evaluate()
    print "jac=",J.getOutput().nz[0]-exp(1)
    self.assertAlmostEqual(J.getOutput()[0,0],exp(1),5,"Evaluation output mismatch")
    
  def test_eval(self):
    self.message('CVodes integration: evaluation')
    num=self.num
    qe=self.qe
    qe.setInput([num['q0']],0)
    qe.setInput([num['p']],1)
    qe.evaluate()

    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(qe.getOutput()[0],q0*exp(tend**3/(3*p)),9,"Evaluation output mismatch")
    
  def test_eval_time_offset(self):
    self.message('CVodes integration: evaluation time offset')
    num=self.num
    integrator=Integrator(self.integrator)
    integrator.setOption("fsens_err_con", True)
    integrator.setOption("t0",0.7)
    integrator.init()

    integrator.setInput([num['q0']],"x0")
    integrator.setInput([num['p']],"p")
    integrator.evaluate()

    tend=num['tend']
    q0=num['q0']
    p=num['p']

    self.assertAlmostEqual(integrator.getOutput()[0],q0*exp((tend**3-0.7**3)/(3*p)),9,"Evaluation output mismatch")
    
    
  def test_jac1(self):
    self.message('CVodes integration: jacobian to q0')
    num=self.num
    J=self.qe.jacobian(0)
    J.init()
    J.setInput([num['q0']],0)
    J.setInput([num['p']],1)
    J.evaluate()
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(J.getOutput()[0],exp(tend**3/(3*p)),9,"Evaluation output mismatch")
    
  def test_jac2(self):
    self.message('CVodes integration: jacobian to p')
    num=self.num
    J=self.qe.jacobian(1)
    J.init()
    J.setInput([num['q0']],0)
    J.setInput([num['p']],1)
    J.evaluate()
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(J.getOutput()[0],-(q0*tend**3*exp(tend**3/(3*p)))/(3*p**2),9,"Evaluation output mismatch")
    
  def test_bug_repeat(self):
    num={'tend':2.3,'q0':[0,7.1,7.1],'p':2}
    self.message("Bug that appears when rhs contains repeats")
    A=array([1,0.1,1])
    p0 = 1.13
    y0=A[0]
    yc0=dy0=A[1]
    te=0.4

    t=SXElement.sym("t")
    q=SX.sym("y",3,1)
    p=SXElement.sym("p")

    dh = p+q[0]**2
    f=SXFunction(daeIn(x=q,p=p,t=t),daeOut(ode=vertcat([dh ,q[0],dh])))
    f.init()
    
    integrator = Integrator("cvodes",f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    #integrator.setOption("verbose",True)
    integrator.setOption("fsens_err_con", True)
    integrator.setOption("steps_per_checkpoint",10000)
    integrator.setOption("t0",0)
    integrator.setOption("tf",te)

    integrator.init()

    q0   = MX.sym("q0",3,1)
    par  = MX.sym("p",1,1)
    qend, = integratorOut(integrator.call(integratorIn(x0=q0,p=par)),"xf")
    qe=MXFunction([q0,par],[qend])
    qe.init()

    #J=self.qe.jacobian(2)
    J=qe.jacobian(0)
    J.init()
    J.setInput(A,0)
    J.setInput(p0,1)
    J.evaluate()
    outA=J.getOutput().toArray()
    f=SXFunction(daeIn(x=q,p=p,t=t),daeOut(ode=vertcat([dh ,q[0],(1+1e-9)*dh])))
    f.init()
    
    integrator = Integrator("cvodes",f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    #integrator.setOption("verbose",True)
    integrator.setOption("steps_per_checkpoint",10000)
    integrator.setOption("fsens_err_con", True)
    integrator.setOption("t0",0)
    integrator.setOption("tf",te)

    integrator.init()

    q0   = MX.sym("q0",3,1)
    par  = MX.sym("p",1,1)
    qend, = integratorOut(integrator.call(integratorIn(x0=q0,p=par)),"xf")
    qe=MXFunction([q0,par],[qend])
    qe.init()

    #J=self.qe.jacobian(2)
    J=qe.jacobian(0)
    J.init()
    J.setInput(A,0)
    J.setInput(p0,1)
    J.evaluate()
    outB=J.output().toArray()
    print outA-outB
    
  def test_hess3(self):
    self.message('CVodes integration: hessian to p: Jacobian of integrator.jacobian')
    num=self.num
    J=self.integrator.jacobian("p","xf")
    J.init()
    H=J.jacobian("p")
    H.init()
    H.setInput([num['q0']],"x0")
    H.setInput([num['p']],"p")
    H.evaluate()
    num=self.num
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(H.getOutput()[0],(q0*tend**6*exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*exp(tend**3/(3*p)))/(3*p**3),9,"Evaluation output mismatch")

  def test_hess4(self):
    self.message('CVodes integration: hessian to p: Jacobian of integrator.jacobian indirect')
    num=self.num
    J=self.integrator.jacobian("p","xf")
    J.init()
    
    q0=MX.sym("q0")
    p=MX.sym("p")
    Ji = MXFunction([q0,p],J.call(integratorIn(x0=q0,p=p)))
    #Ji.setOption("ad_mode","reverse")
    Ji.init()
    H=Ji.jacobian(1)
    H.init()
    H.setInput([num['q0']],0)
    H.setInput([num['p']],1)
    H.evaluate()
    num=self.num
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(H.getOutput()[0],(q0*tend**6*exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*exp(tend**3/(3*p)))/(3*p**3),9,"Evaluation output mismatch")

  def test_hess5(self):
    self.message('CVodes integration: hessian to p in an MX tree')
    num=self.num
    q0=MX.sym("q0")
    p=MX.sym("p")
    qe = MXFunction([q0,p],self.integrator.call(integratorIn(x0=q0,p=p)))
    qe.init()

    JT = MXFunction([q0,p],[qe.jac(1,0)[0].T])
    JT.init()
    JT.setInput([num['q0']],0)
    JT.setInput([num['p']],1)
    JT.evaluate()
    print JT.getOutput()

    H  = JT.jacobian(1)
    H.init()
    H.setInput([num['q0']],0)
    H.setInput([num['p']],1)
    H.evaluate()
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(H.getOutput()[0],(q0*tend**6*exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*exp(tend**3/(3*p)))/(3*p**3),9,"Evaluation output mismatch")
    
  def test_hess6(self):
    self.message('CVodes integration: hessian to p in an MX tree')
    num=self.num
    q0=MX.sym("q0")
    p=MX.sym("p")
    qe = MXFunction([q0,p],self.integrator.call(integratorIn(x0=q0,p=p)))
    qe.init()
    
    H = qe.hessian(1)
    H.init()
    H.setInput([num['q0']],0)
    H.setInput([num['p']],1)
    H.evaluate()
    num=self.num
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(H.getOutput()[0],(q0*tend**6*exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*exp(tend**3/(3*p)))/(3*p**3),9,"Evaluation output mismatch")
     
  def test_glibcbug(self):
    self.message("former glibc error")
    A=array([2.3,4.3,7.6])
    B=array([[1,2.3,4],[-2,1.3,4.7],[-2,6,9]])

    te=0.7
    t=SX.sym("t")
    q=SX.sym("q",3,1)
    p=SX.sym("p",9,1)
    f_in = daeIn(t=t, x=q, p=p)
    f_out = daeOut(ode=mul(c.reshape(p,3,3),q))
    f=SXFunction(f_in,f_out)
    f.init()
    integrator = Integrator("cvodes",f)
    integrator.setOption("fsens_err_con", True)
    integrator.setOption("steps_per_checkpoint",1000)
    integrator.setOption("t0",0)
    integrator.setOption("tf",te)
    integrator.init()
    q0   = MX.sym("q0",3,1)
    par  = MX.sym("p",9,1)
    qend, = integratorOut(integrator.call(integratorIn(x0=q0,p=par)),"xf")
    qe=integrator.jacobian("p","xf")
    qe.init()
    qe = qe.call(integratorIn(x0=q0,p=par))[0]

    qef=MXFunction([q0,par],[qe])
    qef.init()

    qef.setInput(A,0)
    qef.setInput(B.ravel(),1)
    qef.evaluate()
    
  def test_linear_system(self):
    self.message("Linear ODE")
    if not(scipy_available):
        return
    A=array([2.3,4.3,7.6])
    B=array([[1,2.3,4],[-2,1.3,4.7],[-2,6,9]])
    te=0.7
    Be=expm(B*te)
    t=SX.sym("t")
    q=SX.sym("q",3,1)
    p=SX.sym("p",9,1)

    f=SXFunction(daeIn(t=t,x=q,p=p),daeOut(ode=mul(c.reshape(p,3,3),q)))
    f.init()

    integrator = Integrator("cvodes",f)
    integrator.setOption("fsens_err_con", True)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    #integrator.setOption("verbose",True)
    integrator.setOption("steps_per_checkpoint",10000)
    integrator.setOption("t0",0)
    integrator.setOption("tf",te)

    integrator.init()

    q0   = MX.sym("q0",3,1)
    par  = MX.sym("p",9,1)
    qend, = integratorOut(integrator.call(integratorIn(x0=q0,p=par)),"xf")
    qe=MXFunction([q0,par],[qend])
    qe.init()
    qendJ=integrator.jacobian("x0","xf")
    qendJ.init()
    qendJ = qendJ.call(integratorIn(x0=q0,p=par))[0]

    qeJ=MXFunction([q0,par],[qendJ])
    qeJ.init()

    qendJ2=integrator.jacobian("x0","xf")
    qendJ2.init()
    qendJ2 = qendJ2.call(integratorIn(x0=q0,p=par))[0]

    qeJ2=MXFunction([q0,par],[qendJ2])
    qeJ2.init()
    
    qe.setInput(A,0)
    qe.setInput(vec(B),1)
    qe.evaluate()
    self.checkarray(dot(Be,A)/1e3,qe.getOutput()/1e3,"jacobian(INTEGRATOR_X0,INTEGRATOR_XF)")
    qeJ.setInput(A,0)
    qeJ.setInput(vec(B),1)
    qeJ.evaluate()
    self.checkarray(qeJ.getOutput()/1e3,Be/1e3,"jacobian(INTEGRATOR_X0,INTEGRATOR_XF)")
    
    
    qeJ2.setInput(A,0)
    qeJ2.setInput(vec(B),1)
    qeJ2.evaluate()
    
    return # this should return identical zero
    H=qeJ.jacobian(0,0)
    H.setOption("ad_mode","reverse")
    H.init()
    H.setInput(A,0)
    H.setInput(vec(B),1)
    H.evaluate()
    print array(H.getOutput())
    
    
  def test_mathieu_system(self):
    self.message("Mathieu ODE")
    A=array([0.3,1.2])
    B=array([1.3,4.3,2.7])
    te=0.7

    t=SX.sym("t")
    q=SX.sym("q",2,1)
    p=SX.sym("p",3,1)

    f=SXFunction(daeIn(x=q,p=p,t=t),daeOut(ode=vertcat([q[1],(p[0]-2*p[1]*cos(2*p[2]))*q[0]])))
    f.init()
    
    integrator = Integrator("cvodes",f)
    integrator.setOption("fsens_err_con", True)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    #integrator.setOption("verbose",True)
    integrator.setOption("steps_per_checkpoint",10000)
    integrator.setOption("t0",0)
    integrator.setOption("tf",te)

    integrator.init()

    q0   = MX.sym("q0",2,1)
    par  = MX.sym("p",3,1)
    qend, = integratorOut(integrator.call(integratorIn(x0=q0,p=par)),"xf")
    qe=MXFunction([q0,par],[qend])
    qe.init()
    qendJ=integrator.jacobian("x0","xf")
    qendJ.init()
    qendJ =qendJ.call(integratorIn(x0=q0,p=par))[0]
    qeJ=MXFunction([q0,par],[qendJ])
    qeJ.init()

    qe.setInput(A,0)
    qe.setInput(B,1)
    qe.evaluate()
    print array(qe.getOutput())

  def test_nl_system(self):
    """
    y'' = a + (y')^2 , y(0)=y0, y'(0)=yc0
    
    The solution is:
    y=(2*y0-log(yc0^2/a+1))/2-log(cos(atan(yc0/sqrt(a))+sqrt(a)*t))

    """
    self.message("Nonlinear ODE sys")
    A=array([1,0.1])
    p0 = 1.13
    y0=A[0]
    yc0=dy0=A[1]
    te=0.4

    t=SX.sym("t")
    q=SX.sym("y",2,1)
    p=SX.sym("p",1,1)
    # y
    # y'
    f=SXFunction(daeIn(x=q,p=p,t=t),daeOut(ode=vertcat([q[1],p[0]+q[1]**2 ])))
    f.init()
    
    integrator = Integrator("cvodes",f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    #integrator.setOption("verbose",True)
    integrator.setOption("steps_per_checkpoint",10000)
    integrator.setOption("fsens_err_con", True)
    integrator.setOption("t0",0)
    integrator.setOption("tf",te)

    integrator.init()

    t0   = MX(0)
    tend = MX(te)
    q0   = MX.sym("q0",2,1)
    par  = MX.sym("p",1,1)
    qend, = integratorOut(integrator.call(integratorIn(x0=q0,p=par)),"xf")
    qe=MXFunction([q0,par],[qend])
    qe.init()
    qendJ=integrator.jacobian("x0","xf")
    qendJ.init()
    qendJ = qendJ.call(integratorIn(x0=q0,p=par))[0]
    qeJ=MXFunction([q0,par],[qendJ])
    qeJ.init()

    qe.setInput(A,0)
    qe.setInput(p0,1)
    qe.evaluate()

    print qe.getOutput()[0]
    print qe.getOutput()[1]
    
    self.assertAlmostEqual(qe.getOutput()[0],(2*y0-log(yc0**2/p0+1))/2-log(cos(arctan(yc0/sqrt(p0))+sqrt(p0)*te)),11,"Nonlin ODE")
    self.assertAlmostEqual(qe.getOutput()[1],sqrt(p0)*tan(arctan(yc0/sqrt(p0))+sqrt(p0)*te),11,"Nonlin ODE")
    
    qeJ.setInput(A,0)
    qeJ.setInput(p0,1)
    qeJ.evaluate()
    
    Jr = array([[1,(sqrt(p0)*tan(sqrt(p0)*te+arctan(dy0/sqrt(p0)))-dy0)/(dy0**2+p0)],[0,(p0*tan(sqrt(p0)*te+arctan(dy0/sqrt(p0)))**2+p0)/(dy0**2+p0)]])
    self.checkarray(qeJ.getOutput(),Jr,"jacobian of Nonlin ODE")
    
    #qe.setOption("ad_mode","reverse")
    #qe.init()
    Jf=qe.jacobian(0,0)
    Jf.init()
    Jf.setInput(A,0)
    Jf.setInput(p0,1)
    Jf.evaluate()
    self.checkarray(Jf.getOutput(),Jr,"Jacobian of Nonlin ODE")
    
    #qe.setOption("ad_mode","forward")
    #qe.init();
    Jf=qe.jacobian(0,0)
    Jf.init()
    Jf.setInput(A,0)
    Jf.setInput(p0,1)
    Jf.evaluate()
    self.checkarray(Jf.getOutput(),Jr,"Jacobian of Nonlin ODE")
    
    # Joel: This is no longer supported: might be a good idea to support again, though
    #qeJ=integrator.jac("x0","xf")
    #qeJ.init()
    #qeJ.setInput(list(A)+[0,1,0,0],"x0")
    #qeJ.setAdjSeed([0,0]+[0,1,0,0],"xf")
    #qeJ.evaluate(0,1)
    #print qeJ.getOutput()
    #print qeJ.getAdjSens("x0")
    
    Jr = matrix([[(sqrt(p0)*(te*yc0**2-yc0+p0*te)*tan(arctan(yc0/sqrt(p0))+sqrt(p0)*te)+yc0**2)/(2*p0*yc0**2+2*p0**2)],[(sqrt(p0)*((te*yc0**2-yc0+p0*te)*tan(arctan(yc0/sqrt(p0))+sqrt(p0)*te)**2+te*yc0**2-yc0+p0*te)+(yc0**2+p0)*tan(arctan(yc0/sqrt(p0))+sqrt(p0)*te))/(sqrt(p0)*(2*yc0**2+2*p0))]])  
    
    #qe.setOption("ad_mode","reverse")
    #qe.init()
    Jf=qe.jacobian(1,0)
    Jf.init()
    Jf.setInput(A,0)
    Jf.setInput(p0,1)
    Jf.evaluate()
    self.checkarray(Jf.getOutput(),Jr,"Jacobian of Nonlin ODE")
    
    #qe.setOption("ad_mode","forward")
    #qe.init()
    Jf=qe.jacobian(1,0)
    Jf.init()
    Jf.setInput(A,0)
    Jf.setInput(p0,1)
    Jf.evaluate()
    self.checkarray(Jf.getOutput(),Jr,"Jacobian of Nonlin ODE")
    
    qendJ=integrator.jacobian("p","xf")
    qendJ.init()
    qendJ = qendJ.call(integratorIn(x0=q0,p=par))[0]
    qeJ=MXFunction([q0,par],[qendJ])
    qeJ.init()

    qeJ.setInput(A,0)
    qeJ.setInput(p0,1)
    qeJ.evaluate()
    
    self.checkarray(qeJ.getOutput(),Jr,"jacobian of Nonlin ODE")
    
    
    
    
    qeJf=MXFunction([q0,par],[vec(qeJ.call([q0,par])[0])])
    #qeJf.setOption("ad_mode","reverse")
    qeJf.init()
    
    H=qeJf.jacobian(0,0)
    H.init()
    H.setInput(A,0)
    H.setInput(p0,1)
    H.evaluate()
    def sec(x):
      return 1.0/cos(x)
    Hr = array([[0,0],[0,-(2*yc0*tan(arctan(yc0)+te))/(yc0**4+2*yc0**2+1)+sec(arctan(yc0)+te)**2/(yc0**4+2*yc0**2+1)+(2*yc0**2)/(yc0**4+2*yc0**2+1)-1/(yc0**2+1)],[0,0],[0,-(2*yc0*tan(arctan(yc0)+te)**2)/(yc0**4+2*yc0**2+1)+(2*sec(arctan(yc0)+te)**2*tan(arctan(yc0)+te))/(yc0**4+2*yc0**2+1)-(2*yc0)/(yc0**4+2*yc0**2+1)]])
    print array(H.getOutput())
    print Hr
        

  def test_hessian2D(self):
    self.message("hessian")
    N=2

    x0_ = DMatrix([1,0.1])
    A_  = DMatrix([[3,1],[0.74,4]])

    A = SX.sym("A",N,N)
    x = SX.sym("x",N)

    ode = SXFunction(daeIn(x=x, p=vec(A)),daeOut(ode=mul(A,x)))
    I = Integrator("cvodes",ode)
    I.setOption("fsens_err_con", True)
    I.setOption('reltol',1e-12)
    I.init()
    I.setInput(x0_,"x0")
    I.setInput(vec(A_),"p")
    I.evaluate()

    q0=MX.sym("q0",N)
    p=MX.sym("p",N*N)
    qe = MXFunction([q0,p],I.call(integratorIn(x0=q0,p=p)))
    qe.init()

    JT = MXFunction([q0,p],[qe.jac(1,0).T])
    JT.init()

    H  = JT.jacobian(1)
    H.init()
    H.setInput(x0_,0)
    H.setInput(vec(A_),1)
    H.evaluate()

    H1 = H.getOutput()
    
    ## Joel: Only Hessians of scalar functions allowed
    #H = qe.hessian(1)
    #H.init()
    #H.setInput(x0_,0)
    #H.setInput(vec(A_),1)
    #H.evaluate()
    #H2 = H.getOutput()
    
    #self.checkarray(H1,H2,"hessian")
    
  def test_issue535(self):
    self.message("regression test for #535")
    t=SX.sym("t")
    x=SX.sym("x")
    rx=SX.sym("rx")
    p=SX.sym("p")
    dp=SX.sym("dp")

    z=SX.sym("z")
    rz=SX.sym("rz")
    rp=SX.sym("rp")
    f = SXFunction(daeIn(**{'x': x, 'z': z}),daeOut(**{'alg': x-z, 'ode': z}))
    f.init()
    g = SXFunction(rdaeIn(**{'x': x, 'z': z, 'rx': rx, 'rz': rz}),rdaeOut(**{'alg': x-rz, 'ode': rz}))
    g.init()

    integrator = Integrator("idas",f,g)
    integrator.setOption({'calc_ic': True, 'tf': 2.3, 'reltol': 1e-10, 'augmented_options': {'reltol': 1e-09, 'abstol': 1e-09 }, 'calc_icB': True, 'abstol': 1e-10, 't0': 0.2})
    integrator.init()

    integrator.setInput(7.1,"x0")
    if not integrator.getInput("p").isEmpty():
      integrator.setInput(2,"p")
    if not integrator.getInput("rx0").isEmpty():
      integrator.setInput(0.13,"rx0")
    if not integrator.getInput("rp").isEmpty():
      integrator.setInput(0.127,"rp")

    integrator.evaluate()
    
  def test_collocationPoints(self):
    self.message("collocation points")
    with self.assertRaises(Exception):
      collocationPoints(0,"radau")
    with self.assertRaises(Exception): 
      collocationPoints(10,"radau")
    with self.assertRaises(Exception):
      collocationPoints(0,"legendre")
    with self.assertRaises(Exception): 
      collocationPoints(10,"legendre")
    with self.assertRaises(Exception):
      collocationPoints(1,"foo")
      
    for k in range(1,10):
      r = collocationPoints(k,"radau")
      self.assertEqual(len(r),k+1)
      self.checkarray(DMatrix(r[-1]),DMatrix([1]))
    for k in range(1,10):
      r = collocationPoints(k,"legendre")
      self.assertEqual(len(r),k+1) 
      
if __name__ == '__main__':
    unittest.main()

