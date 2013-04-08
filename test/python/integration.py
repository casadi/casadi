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
  integrators.append((CVodesIntegrator,["ode"],{"abstol": 1e-15,"reltol":1e-15,"fsens_err_con": True,"quad_err_con": False}))
except:
  pass
  
try:
  integrators.append((IdasIntegrator,["dae","ode"],{"abstol": 1e-15,"reltol":1e-15,"fsens_err_con": True,"calc_icB":True}))
except:
  pass

integrators.append((CollocationIntegrator,["dae","ode"],{"implicit_solver":KinsolSolver,"number_of_finite_elements": 18,"startup_integrator":CVodesIntegrator}))
#integrators.append((CollocationIntegrator,["dae","ode"],{"implicit_solver":NLPImplicitSolver,"number_of_finite_elements": 100,"startup_integrator":CVodesIntegrator,"implicit_solver_options": {"nlp_solver": IpoptSolver,"linear_solver_creator": CSparse}}))
#integrators.append((RKIntegrator,["ode"],{"number_of_finite_elements": 1000}))

print "Will test these integrators:"
for cl, t, options in integrators:
  print cl.__name__, " : ", t
  
class Integrationtests(casadiTestCase):

  @skip(memcheck)
  def test_jac(self):
    self.message("Test exact jacobian #536")
    # This test is not automized, but works by inspection only.
    # To activate, recompile after ucnommenting the printout lines in cvodes.c, near "Used for validating casadi#536"
    #return
    DMatrix.setPrecision(18)

    tstart = ssym("tstart")
    tend = ssym("tend")
    
    integrators = [
              (IdasIntegrator,["dae","ode"],{"abstol": 1e-9,"reltol":1e-9,"fsens_err_con": True,"calc_ic":True,"calc_icB":True}),
              (CVodesIntegrator,["ode"],{"abstol": 1e-5,"reltol":1e-5,"fsens_err_con": False,"quad_err_con": False})
              ]

    def variations(p_features, din, dout, rdin, rdout, *args):
      if "ode" in p_features:
        p_features_ = copy.copy(p_features)
        p_features_[p_features.index("ode")] = "dae"
        din_ = copy.copy(din)
        dout_ = copy.copy(dout)
        rdin_ = copy.copy(rdin)
        rdout_ = copy.copy(rdout)
        z = ssym("x", din_["x"].shape)
        din_["z"] = z
        dout_["ode"] = z
        dout_["alg"] = ( dout["ode"] - z) * (-0.8)
        if len(rdin_)>0:
          rz = ssym("rx", rdin_["rx"].shape)
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
      
      x  = ssym("x")
      rx = ssym("rx")
      t = ssym("t")

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
          self.message(Integrator.__name__)
          if p_features[0] in features:
            g = FX()
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
              yield {"iterative_solver"+post: "tfqmr", "use_preconditionerB": True, "linear_solverB" : CSparse} # Bug in Sundials? Preconditioning seems to be needed
             
            def solveroptions(post=""):
              yield {"linear_solver_type" +post: "dense" }
              #for it in itoptions(post):
              #  d = {"linear_solver_type" +post: "iterative" }
              #  d.update(it)
              #  yield d
              #yield {"linear_solver_type" +post: "banded", "lower_bandwidth"+post: 0, "upper_bandwidth"+post: 0 }
              yield {"linear_solver_type" +post: "user_defined", "linear_solver"+post: CSparse }
                
            for a_options in solveroptions("B"):
              for f_options in solveroptions():
                message = "f_options: %s , a_options: %s" % (str(f_options) , str(a_options))
                print message
                integrator = Integrator(f,g)
                integrator.setOption("exact_jacobianB",True)
                integrator.setOption("gather_stats",True)
                integrator.setOption("verbose",True)
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
                    if not ff.input(i).empty():
                      ff.input(i).set(v)

                integrator.evaluate(0,0)
                fs.evaluate(0,0)
                print "res=",integrator.output("xf")-fs.output("xf"), fs.output("xf")
                print "Rres=",integrator.output("rxf")-fs.output("rxf"), fs.output("rxf")
                # self.checkarray(integrator.output("rxf"),fs.output("rxf"),digits=4)
                stats = integrator.getStats()
                
                print stats
                self.assertTrue(stats["nsteps"]<1500)
                self.assertTrue(stats["nstepsB"]<2500)
                self.assertTrue(stats["nlinsetups"]<100)
                self.assertTrue(stats["nlinsetupsB"]<250)

  @skip(memcheck)
  def test_lsolvers(self):
    self.message("Test different linear solvers")

    tstart = ssym("tstart")
    tend = ssym("tend")
    
    integrators = [
              (IdasIntegrator,["dae","ode"],{"abstol": 1e-9,"reltol":1e-9,"fsens_err_con": True,"calc_ic":True,"calc_icB":True}),
              (CVodesIntegrator,["ode"],{"abstol": 1e-15,"reltol":1e-15,"fsens_err_con": True,"quad_err_con": False})
              ]
              
    def checks():  
      t=SX("t")
      x=SX("x")
      rx=SX("rx")
      p=SX("p")
      dp=SX("dp")

      z=SX("z")
      rz=SX("rz")
      rp=SX("rp")    
      solutionin = {'x0':x, 'p': p, 'rx0': rx,'rp' : rp}            
      pointA = {'x0':7.1,'p': 2, 'rx0': 0.13, 'rp': 0.127}
      ti = (0.2,2.3)
      yield (["dae"],{'x': x, 'z': z},{'alg': x-z, 'ode': z},{'x': x, 'z': z, 'rx': rx, 'rz': rz},{'alg': x-rz, 'ode': rz},solutionin,{'rxf': rx+x*(exp(tend-tstart)-1), 'xf':x*exp(tend-tstart)},pointA,ti)
      yield (["dae"],{'x': x, 'z': z},{'alg': x-z, 'ode': z},{'x': x, 'z': z, 'rx': rx, 'rz': rz},{'alg': rx-rz, 'ode': rz},solutionin,{'rxf': rx*exp(tend-tstart), 'xf':x*exp(tend-tstart)},pointA,ti)
      yield (["ode"],{'x': x},{'ode': x},{'x': x,'rx': rx},{'ode': x},solutionin,{'rxf': rx+x*(exp(tend-tstart)-1), 'xf':x*exp(tend-tstart)},pointA,ti)
      yield (["ode"],{'x': x},{'ode': x},{'x': x,'rx': rx},{'ode': rx},solutionin,{'rxf': rx*exp(tend-tstart), 'xf':x*exp(tend-tstart)},pointA,ti)
      
      A=array([1,0.1])
      p0 = 1.13

      q=ssym("y",2,1)
      y0=q[0]
      yc0=dy0=q[1]
      p=ssym("p",1,1)

      s1=(2*y0-log(yc0**2/p+1))/2-log(cos(arctan(yc0/sqrt(p))+sqrt(p)*(tend-tstart)))
      s2=sqrt(p)*tan(arctan(yc0/sqrt(p))+sqrt(p)*(tend-tstart))
      yield (["ode"],{'x':q,'p':p},{'ode': vertcat([q[1],p[0]+q[1]**2 ])},{},{},{'x0':q, 'p': p} ,{'xf': vertcat([s1,s2])},{'x0': A, 'p': p0},(0,0.4) )

    for p_features, din, dout, rdin, rdout, solutionin, solution, point, (tstart_, tend_) in checks():

      for Integrator, features, options in integrators:
        self.message(Integrator.__name__)
        if p_features[0] in features:
          g = FX()
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
            yield {"iterative_solver"+post: "tfqmr", "use_preconditionerB": True, "linear_solverB" : CSparse} # Bug in Sundials? Preconditioning seems to be needed
           
          def solveroptions(post=""):
            yield {"linear_solver_type" +post: "dense" }
            for it in itoptions(post):
              d = {"linear_solver_type" +post: "iterative" }
              d.update(it)
              yield d
            #yield {"linear_solver_type" +post: "banded", "lower_bandwidth"+post: 0, "upper_bandwidth"+post: 0 }
            yield {"linear_solver_type" +post: "user_defined", "linear_solver"+post: CSparse }
              
          for a_options in solveroptions("B"):
            for f_options in solveroptions():
              message = "f_options: %s , a_options: %s" % (str(f_options) , str(a_options))
              print message
              integrator = Integrator(f,g)
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
                  if not ff.input(i).empty():
                    ff.input(i).set(v)

              integrator.evaluate(1,0)
              
              self.checkfx(integrator,fs,gradient=False,hessian=False,sens_der=False,digits=4,digits_sens=4,failmessage=message,verbose=False)
              
              


  @skip(memcheck)
  def test_X(self):
    self.message("Extensive integrator tests")
    
    num=self.num
    tstart = ssym("tstart")
    tend = ssym("tstart")

    
    for Integrator, features, options in integrators:
      self.message(Integrator.__name__)
        
        
      def variations(p_features, din, dout, rdin, rdout, *args):
        if "ode" in p_features:
          p_features_ = copy.copy(p_features)
          p_features_[p_features.index("ode")] = "dae"
          din_ = copy.copy(din)
          dout_ = copy.copy(dout)
          rdin_ = copy.copy(rdin)
          rdout_ = copy.copy(rdout)
          z = ssym("x", din_["x"].shape)
          din_["z"] = z
          dout_["ode"] = z
          dout_["alg"] = ( dout["ode"] - z) * (-0.8)
          if len(rdin_)>0:
            rz = ssym("rx", rdin_["rx"].shape)
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
        t=ssym("t")
        x=ssym("x")
        rx=ssym("rx")
        p=ssym("p")
        dp=ssym("dp")

        z=ssym("z")
        rz=ssym("rz")
        rp=ssym("rp")
        
        si = {'x0':x, 'p': p, 'rx0': rx,'rp' : rp}            
        pointA = {'x0':x0,'p': p_, 'rx0': rx0_, 'rp': 0.127}
        
        ti = (0.2,num['tend'])
        yield (["ode"],{'x':x},{'ode': 0},{},{},si,{'xf':x},pointA,ti)
        yield (["ode"],{'x':x},{'ode': 1},{},{},si,{'xf':x+(tend-tstart)},pointA,ti)
        yield (["ode"],{'x':x},{'ode': x},{},{},si,{'xf':x*exp(tend-tstart)},pointA,ti)
        yield (["ode"],{'x':x,'t':t},{'ode': t},{},{},si,{'xf':x+(tend**2/2-tstart**2/2)},pointA,ti)
        yield (["ode"],{'x':x,'t':t},{'ode': x*t},{},{},si,{'xf':x*exp(tend**2/2-tstart**2/2)},pointA,ti)
        yield (["ode"],{'x':x,'p':p},{'ode': x/p},{},{},si,{'xf':x*exp((tend-tstart)/p)},pointA,ti)
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

        q=ssym("y",2,1)
        y0=q[0]
        yc0=dy0=q[1]
        p=ssym("p",1,1)
        
        s1=(2*y0-log(yc0**2/p+1))/2-log(cos(arctan(yc0/sqrt(p))+sqrt(p)*(tend-tstart)))
        s2=sqrt(p)*tan(arctan(yc0/sqrt(p))+sqrt(p)*(tend-tstart))
        yield (["ode"],{'x':q,'p':p},{'ode': vertcat([q[1],p[0]+q[1]**2 ])},{},{},{'x0':q, 'p': p} ,{'xf': vertcat([s1,s2])},{'x0': A, 'p': p0},(0,0.4) )
      
      for tt in checks():
        print tt
        for p_features, din, dout, rdin, rdout, solutionin, solution, point, (tstart_, tend_) in variations(*tt):
          if p_features[0] in features:
            message = "%s: %s => %s, %s => %s, explicit (%s) tstart = %f" % (Integrator.__name__,str(din),str(dout),str(rdin),str(rdout),str(solution),tstart_)
            print message
            g = FX()
            if len(rdin)>1:
              g = SXFunction(rdaeIn(**rdin),rdaeOut(**rdout))
              g.init()
               
            f = SXFunction(daeIn(**din),daeOut(**dout))
            f.init()
            
            for k in solution.keys():
              solution[k] = substitute(solution[k],vertcat([tstart,tend]),vertcat([tstart_,tend_]))
            
            fs = SXFunction(integratorIn(**solutionin),integratorOut(**solution))
            fs.init()
            
            integrator = Integrator(f,g)
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
#t=SX("t")
#x=SX("x")
#rx=SX("rx")
#p=SX("p")
#dp=SX("dp")

#z=SX("z")
#rz=SX("rz")
#rp=SX("rp")
#f = SXFunction(daeIn(**{din}),daeOut(**{dout}))
#f.init()
#g = SXFunction(rdaeIn(**{rdin}),rdaeOut(**{rdout}))
#g.init()

#integrator = {intclass.__name__}(f,g)
#integrator.setOption({options})
#integrator.init()

#integrator.input("x0").set({x0})
#if not integrator.input("p").empty():
#  integrator.input("p").set({p_})
#if not integrator.input("rx0").empty():
#  integrator.input("rx0").set(0.13)
#if not integrator.input("rp").empty():
#  integrator.input("rp").set(0.127)
#              """.format(din=din,dout=dout,rdin=rdin,rdout=rdout,x0=x0,p_=p_,intclass=Integrator,options=integrator.dictionary())
#              message+="\nTo reproduce:\n" + reproduce

                
       
            for ff in [fs,integrator]:
              for k,v in point.items():
                i = getattr(casadi,('integrator_'+k).upper())
                if not ff.input(i).empty():
                  ff.input(i).set(v)
            integrator.evaluate(1,0)
            
            self.checkfx(integrator,fs,gradient=False,hessian=False,sens_der=False,digits=4,digits_sens=4,failmessage=message,verbose=False)

        
  def setUp(self):
    # Reference solution is x0 e^((t^3-t0^3)/(3 p))
    t=ssym("t")
    x=ssym("x")
    p=ssym("p")
    f=SXFunction(daeIn(t=t, x=x, p=p),daeOut(ode=x/p*t**2))
    f.init()
    integrator = CVodesIntegrator(f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    integrator.setOption("fsens_err_con", True)
    #integrator.setOption("verbose",True)
    integrator.setOption("t0",0)
    integrator.setOption("tf",2.3)
    integrator.init()
    q0   = MX("q0")
    par  = MX("p")
    
    # qend,*_ = integrator.call([q0,par]) # Valid Python3 syntax
    qend,_,_,_ = integrator.call([q0,par])
    
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
  
    
  def test_parameterize_time(self):
    self.message("parametrizeTime")
    num=self.num
    f = self.f
    
    f_ = parameterizeTime(f)
    for intf in [CVodesIntegrator, IdasIntegrator]:
      integrator = intf(parameterizeTime(f))
      integrator.setOption("reltol",1e-15)
      integrator.setOption("abstol",1e-15)
      integrator.setOption("fsens_err_con", True)
      #integrator.setOption("verbose",True)
      integrator.setOption("t0",0)
      integrator.setOption("tf",1)
      integrator.init()
      
      tend = num['tend']
      p = num['p']
      q0 = num['q0']
      
      t0 = 0.7
      
      integrator.input("p").set([0,tend,p])
      integrator.input("x0").set([q0])
      
      integrator.evaluate()
      
      self.assertAlmostEqual(integrator.output()[0],q0*exp(tend**3/(3*p)),9,"Evaluation output mismatch")
    
      # Integrate with time offset
      integrator.input("p").set([t0,tend,p])
      integrator.input("x0").set([q0])
      
      integrator.evaluate()
      
      self.assertAlmostEqual(integrator.output()[0],q0*exp((tend**3-t0**3)/(3*p)),9,"Evaluation output mismatch")
      
      # Forward sensitivity to q0
      integrator.input("x0").set([q0])
      integrator.fwdSeed("x0").set(1)
      integrator.evaluate(1,0)
      
      self.assertAlmostEqual(integrator.fwdSens()[0],exp((tend**3-t0**3)/(3*p)),9,"Evaluation output mismatch")
      
      # Forward sensitivity to p
      integrator.fwdSeed("x0").set(0)
      integrator.fwdSeed("p").set([0,0,1])
      integrator.evaluate(1,0)
      
      self.assertAlmostEqual(integrator.fwdSens()[0],-(q0*(tend**3-t0**3)*exp((tend**3-t0**3)/(3*p)))/(3*p**2),9,"Evaluation output mismatch")
      
      # Forward sensitivity to tf
      integrator.fwdSeed("x0").set(0)
      integrator.fwdSeed("p").set([0,1,0])
      integrator.evaluate(1,0)
      
      self.assertAlmostEqual(integrator.fwdSens()[0],(q0*tend**2*exp((tend**3-t0**3)/(3*p)))/p,7,"Evaluation output mismatch")
      
      # Forward sensitivity to t0
      integrator.fwdSeed("x0").set(0)
      integrator.fwdSeed("p").set([1,0,0])
      integrator.input("p").set([t0,tend,p])
      integrator.evaluate(1,0)
      
      self.assertAlmostEqual(integrator.fwdSens()[0],-(q0*t0**2*exp((tend**3-t0**3)/(3*p)))/p,7,"Evaluation output mismatch")

      if not(intf is IdasIntegrator):
        # (*) IDAS backward sens seems to fail for somewhat small tolerances
        integrator.adjSeed("xf").set(1)
        integrator.input("p").set([t0,tend,p])
        integrator.evaluate(0,1)

        self.assertAlmostEqual(integrator.adjSens("x0")[0],exp((tend**3-t0**3)/(3*p)),9,"Evaluation output mismatch")
        self.assertAlmostEqual(integrator.adjSens("p")[2],-(q0*(tend**3-t0**3)*exp((tend**3-t0**3)/(3*p)))/(3*p**2),9,"Evaluation output mismatch")
        self.assertAlmostEqual(integrator.adjSens("p")[1],(q0*tend**2*exp((tend**3-t0**3)/(3*p)))/p,7,"Evaluation output mismatch")
        self.assertAlmostEqual(integrator.adjSens("p")[0],-(q0*t0**2*exp((tend**3-t0**3)/(3*p)))/p,7,"Evaluation output mismatch")
    
      # (*) Try IDAS again with very low tolerances
      if 0:
        integrator = IdasIntegrator(f_)
        integrator.setOption("reltol",1e-6)
        integrator.setOption("abstol",1e-6)
        integrator.setOption("fsens_err_con", True)
        integrator.setOption("t0",0)
        integrator.setOption("tf",1)
        integrator.init()
        
        integrator.adjSeed("x0").set(1)
        integrator.input("x0").set([q0])
        integrator.input("p").set([t0,tend,p])
        integrator.evaluate(0,1)
        self.assertAlmostEqual(integrator.adjSens("x0")[0],exp((tend**3-t0**3)/(3*p)),2,"Evaluation output mismatch")
        print integrator.adjSens("p")[2],-(q0*(tend**3-t0**3)*exp((tend**3-t0**3)/(3*p)))/(3*p**2)
        self.assertAlmostEqual(integrator.adjSens("p")[2],-(q0*(tend**3-t0**3)*exp((tend**3-t0**3)/(3*p)))/(3*p**2),2,"Evaluation output mismatch")
        self.assertAlmostEqual(integrator.adjSens("p")[1],(q0*tend**2*exp((tend**3-t0**3)/(3*p)))/p,2,"Evaluation output mismatch")
        self.assertAlmostEqual(integrator.adjSens("p")[0],-(q0*t0**2*exp((tend**3-t0**3)/(3*p)))/p,2,"Evaluation output mismatch")
      
    
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
    f.input(0).set([num['q0']])
    f.evaluate()

    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(f.output()[0],q0*exp(tend**3/(3*p)),9,"Evaluation output mismatch")
  
  def test_issue92c(self):
    self.message("regression check for issue 92")
    t=SX("t")
    x=SX("x")
    y=SX("y")
    z=x*exp(t)
    f=SXFunction(daeIn(t=t, x=vertcat([x,y])),[vertcat([z,z])])
    f.init()
    # Pass inputs
    f.setInput(1.0,"t")
    f.setInput([1.0,0.0],"x")
    # Pass adjoint seeds
    f.setAdjSeed([1.0,0.0])
    # Evaluate with adjoint mode AD
    f.evaluate(0,1)
    # print result
    print f.output()
    print f.adjSens("x")
  
  def test_issue92b(self):
    self.message("regression check for issue 92")
    t=SX("t")
    x=SX("x")
    y=SX("y")
    f=SXFunction(daeIn(t=t, x=vertcat([x,y])),daeOut(ode=vertcat([x,(1+1e-9)*x])))
    integrator = CVodesIntegrator(f)
    integrator.setOption("fsens_err_con", True)
    integrator.setOption("t0",0)
    integrator.setOption("tf",1)
    integrator.init()
    # Pass inputs
    integrator.setInput([1,0],"x0")
    # Pass adjoint seeds
    integrator.setAdjSeed([1.0,0.0],"xf")
    ## Integrate and calculate sensitivities
    integrator.evaluate(0,1)
    # print result
    print integrator.output("xf")
    print integrator.adjSens("x0")
    
  def test_issue92(self):
    self.message("regression check for issue 92")
    t=SX("t")
    x=SX("x")
    var = MX("var",2,1)

    q = vertcat([x,SX("problem")])

    dq=vertcat([x,x])
    f=SXFunction(daeIn(t=t,x=q),daeOut(ode=dq))
    f.init()

    integrator = CVodesIntegrator(f)
    integrator.setOption("fsens_err_con", True)
    integrator.setOption("reltol",1e-12)
    integrator.setOption("t0",0)
    integrator.setOption("tf",1)
    integrator.init()

    qend,_,_,_ = integrator.call([var])

    f = MXFunction([var],[qend[0]])
    f.init()

    J=f.jacobian(0)
    J.init()
    J.input().set([1,0])
    J.evaluate()
    print "jac=",J.output()[0]-exp(1)
    self.assertAlmostEqual(J.output()[0],exp(1),5,"Evaluation output mismatch")
    
  def test_eval(self):
    self.message('CVodes integration: evaluation')
    num=self.num
    qe=self.qe
    qe.input(0).set([num['q0']])
    qe.input(1).set([num['p']])
    qe.evaluate()

    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(qe.output()[0],q0*exp(tend**3/(3*p)),9,"Evaluation output mismatch")
    
  def test_eval_time_offset(self):
    self.message('CVodes integration: evaluation time offset')
    num=self.num
    integrator=Integrator(self.integrator)
    integrator.setOption("fsens_err_con", True)
    integrator.setOption("t0",0.7)
    integrator.init()

    integrator.input("x0").set([num['q0']])
    integrator.input("p").set([num['p']])
    integrator.evaluate()

    tend=num['tend']
    q0=num['q0']
    p=num['p']

    self.assertAlmostEqual(integrator.output()[0],q0*exp((tend**3-0.7**3)/(3*p)),9,"Evaluation output mismatch")
    
    
  def test_jac1(self):
    self.message('CVodes integration: jacobian to q0')
    num=self.num
    J=self.qe.jacobian(0)
    J.init()
    J.input(0).set([num['q0']])
    J.input(1).set([num['p']])
    J.evaluate()
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(J.output()[0],exp(tend**3/(3*p)),9,"Evaluation output mismatch")
    
  def test_jac2(self):
    self.message('CVodes integration: jacobian to p')
    num=self.num
    J=self.qe.jacobian(1)
    J.init()
    J.input(0).set([num['q0']])
    J.input(1).set([num['p']])
    J.evaluate()
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(J.output()[0],-(q0*tend**3*exp(tend**3/(3*p)))/(3*p**2),9,"Evaluation output mismatch")
    
  def test_bug_repeat(self):
    num={'tend':2.3,'q0':[0,7.1,7.1],'p':2}
    self.message("Bug that appears when rhs contains repeats")
    A=array([1,0.1,1])
    p0 = 1.13
    y0=A[0]
    yc0=dy0=A[1]
    te=0.4

    t=SX("t")
    q=ssym("y",3,1)
    p=SX("p")

    dh = p+q[0]**2
    f=SXFunction(daeIn(x=q,p=p,t=t),daeOut(ode=vertcat([dh ,q[0],dh])))
    f.init()
    
    integrator = CVodesIntegrator(f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    #integrator.setOption("verbose",True)
    integrator.setOption("fsens_err_con", True)
    integrator.setOption("steps_per_checkpoint",10000)
    integrator.setOption("t0",0)
    integrator.setOption("tf",te)

    integrator.init()

    q0   = MX("q0",3,1)
    par  = MX("p",1,1)
    qend,_,_,_ = integrator.call([q0,par])
    qe=MXFunction([q0,par],[qend])
    qe.init()

    #J=self.qe.jacobian(2)
    J=qe.jacobian(0)
    J.init()
    J.input(0).set(A)
    J.input(1).set(p0)
    J.evaluate()
    outA=J.output().toArray()
    f=SXFunction(daeIn(x=q,p=p,t=t),daeOut(ode=vertcat([dh ,q[0],(1+1e-9)*dh])))
    f.init()
    
    integrator = CVodesIntegrator(f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    #integrator.setOption("verbose",True)
    integrator.setOption("steps_per_checkpoint",10000)
    integrator.setOption("fsens_err_con", True)
    integrator.setOption("t0",0)
    integrator.setOption("tf",te)

    integrator.init()

    q0   = MX("q0",3,1)
    par  = MX("p",1,1)
    qend,_,_,_ = integrator.call([q0,par])
    qe=MXFunction([q0,par],[qend])
    qe.init()

    #J=self.qe.jacobian(2)
    J=qe.jacobian(0)
    J.init()
    J.input(0).set(A)
    J.input(1).set(p0)
    J.evaluate()
    outB=J.output().toArray()
    print outA-outB
    
  def test_hess(self):
    self.message('CVodes integration: hessian to p: fwd-over-adjoint on integrator')
    num=self.num
    J=self.integrator.jacobian("p","xf")
    J.setOption("number_of_fwd_dir",0)
    J.setOption("number_of_adj_dir",1)
    J.init()
    J.input("x0").set([num['q0']])
    J.input("p").set([num['p']])
    J.adjSeed().set([1])
    # Evaluate
    J.evaluate(0,1)
      
    tend=num['tend']
    q0=num['q0']
    p=num['p']

    self.assertAlmostEqual(J.adjSens("p")[0],(q0*tend**6*exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*exp(tend**3/(3*p)))/(3*p**3),7,"Evaluation output mismatch")

  def test_hess3(self):
    self.message('CVodes integration: hessian to p: Jacobian of integrator.jacobian')
    num=self.num
    J=self.integrator.jacobian("p","xf")
    J.init()
    H=J.jacobian("p")
    H.init()
    H.input("x0").set([num['q0']])
    H.input("p").set([num['p']])
    H.evaluate(0,0)
    num=self.num
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(H.output()[0],(q0*tend**6*exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*exp(tend**3/(3*p)))/(3*p**3),9,"Evaluation output mismatch")

  def test_hess4(self):
    self.message('CVodes integration: hessian to p: Jacobian of integrator.jacobian indirect')
    num=self.num
    J=self.integrator.jacobian("p","xf")
    J.init()
    
    q0=MX("q0")
    p=MX("p")
    Ji = MXFunction([q0,p],J.call([q0,p]))
    #Ji.setOption("ad_mode","reverse")
    Ji.init()
    H=Ji.jacobian(1)
    H.init()
    H.input(0).set([num['q0']])
    H.input(1).set([num['p']])
    H.evaluate(0,0)
    num=self.num
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(H.output()[0],(q0*tend**6*exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*exp(tend**3/(3*p)))/(3*p**3),9,"Evaluation output mismatch")

  def test_hess5(self):
    self.message('CVodes integration: hessian to p in an MX tree')
    num=self.num
    q0=MX("q0")
    p=MX("p")
    qe = MXFunction([q0,p],self.integrator.call([q0,p]))
    qe.init()

    JT = MXFunction([q0,p],[qe.jac(1,0)[0].T])
    JT.init()
    JT.input(0).set([num['q0']])
    JT.input(1).set([num['p']])
    JT.evaluate(1,0)
    print JT.output()

    H  = JT.jacobian(1)
    H.init()
    H.input(0).set([num['q0']])
    H.input(1).set([num['p']])
    H.evaluate()
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(H.output()[0],(q0*tend**6*exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*exp(tend**3/(3*p)))/(3*p**3),9,"Evaluation output mismatch")
    
  def test_hess6(self):
    self.message('CVodes integration: hessian to p in an MX tree')
    num=self.num
    q0=MX("q0")
    p=MX("p")
    qe = MXFunction([q0,p],self.integrator.call([q0,p]))
    qe.init()
    
    H = qe.hessian(1)
    H.init()
    H.input(0).set([num['q0']])
    H.input(1).set([num['p']])
    H.evaluate()
    num=self.num
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(H.output()[0],(q0*tend**6*exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*exp(tend**3/(3*p)))/(3*p**3),9,"Evaluation output mismatch")
 
  def test_issue87(self):
    return # see issue 87
    self.message('CVodes integration: hessian to p: fwd-over-adjoint on integrator')
    num=self.num
    J=self.qe.jacobian(1)
    J.init()
    J.input(0).set([num['q0']])
    J.input(1).set([num['p']])
    J.fwdSeed(0).set([1])
    J.fwdSeed(1).set([1])
    # Evaluate
    J.evaluate(1,1)
      
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    print (q0*tend**6*exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*exp(tend**3/(3*p)))/(3*p**3)
    print J.adjSens()
    print J.fwdSens()
    self.assertAlmostEqual(J.adjSens(1)[0],(q0*tend**6*exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*exp(tend**3/(3*p)))/(3*p**3),9,"Evaluation output mismatch")
    
    
  def test_glibcbug(self):
    self.message("former glibc error")
    A=array([2.3,4.3,7.6])
    B=array([[1,2.3,4],[-2,1.3,4.7],[-2,6,9]])

    te=0.7
    t=ssym("t")
    q=ssym("q",3,1)
    p=ssym("p",9,1)
    f_in = daeIn(t=t, x=q, p=p)
    f_out = daeOut(ode=mul(c.reshape(p,3,3),q))
    f=SXFunction(f_in,f_out)
    f.init()
    integrator = CVodesIntegrator(f)
    integrator.setOption("fsens_err_con", True)
    integrator.setOption("steps_per_checkpoint",1000)
    integrator.setOption("t0",0)
    integrator.setOption("tf",te)
    integrator.init()
    q0   = MX("q0",3,1)
    par  = MX("p",9,1)
    qend,_,_,_ = integrator.call([q0,par])
    qe=integrator.jacobian("p","xf")
    qe.init()
    qe = qe.call([q0,par])[0]

    qef=MXFunction([q0,par],[qe])
    qef.init()

    qef.input(0).set(A)
    qef.input(1).set(B.ravel())
    qef.evaluate()
    
  def test_linear_system(self):
    self.message("Linear ODE")
    if not(scipy_available):
        return
    A=array([2.3,4.3,7.6])
    B=array([[1,2.3,4],[-2,1.3,4.7],[-2,6,9]])
    te=0.7
    Be=expm(B*te)
    t=ssym("t")
    q=ssym("q",3,1)
    p=ssym("p",9,1)

    f=SXFunction(daeIn(t=t,x=q,p=p),daeOut(ode=mul(c.reshape(p,3,3),q)))
    f.init()

    integrator = CVodesIntegrator(f)
    integrator.setOption("fsens_err_con", True)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    #integrator.setOption("verbose",True)
    integrator.setOption("steps_per_checkpoint",10000)
    integrator.setOption("t0",0)
    integrator.setOption("tf",te)

    integrator.init()

    q0   = MX("q0",3,1)
    par  = MX("p",9,1)
    qend,_,_,_ = integrator.call([q0,par])
    qe=MXFunction([q0,par],[qend])
    qe.init()
    qendJ=integrator.jacobian("x0","xf")
    qendJ.init()
    qendJ = qendJ.call([q0,par])[0]

    qeJ=MXFunction([q0,par],[qendJ])
    qeJ.init()

    qendJ2=integrator.jacobian("x0","xf")
    qendJ2.init()
    qendJ2 = qendJ2.call([q0,par])[0]

    qeJ2=MXFunction([q0,par],[qendJ2])
    qeJ2.init()
    
    qe.input(0).set(A)
    qe.input(1).set(B.ravel())
    qe.evaluate()
    self.checkarray(dot(Be,A)/1e3,qe.output()/1e3,"jacobian(INTEGRATOR_X0,INTEGRATOR_XF)")
    qeJ.input(0).set(A)
    qeJ.input(1).set(B.ravel())
    qeJ.evaluate()
    self.checkarray(qeJ.output()/1e3,Be/1e3,"jacobian(INTEGRATOR_X0,INTEGRATOR_XF)")
    
    
    qeJ2.input(0).set(A)
    qeJ2.input(1).set(B.ravel())
    qeJ2.evaluate()
    
    return # this should return identical zero
    H=qeJ.jacobian(0,0)
    H.setOption("ad_mode","reverse")
    H.init()
    H.input(0).set(A)
    H.input(1).set(B.ravel())
    H.evaluate()
    print array(H.output())
    
    
  def test_mathieu_system(self):
    self.message("Mathieu ODE")
    A=array([0.3,1.2])
    B=array([1.3,4.3,2.7])
    te=0.7

    t=ssym("t")
    q=ssym("q",2,1)
    p=ssym("p",3,1)

    f=SXFunction(daeIn(x=q,p=p,t=t),daeOut(ode=vertcat([q[1],(p[0]-2*p[1]*cos(2*p[2]))*q[0]])))
    f.init()
    
    integrator = CVodesIntegrator(f)
    integrator.setOption("fsens_err_con", True)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    #integrator.setOption("verbose",True)
    integrator.setOption("steps_per_checkpoint",10000)
    integrator.setOption("t0",0)
    integrator.setOption("tf",te)

    integrator.init()

    q0   = MX("q0",2,1)
    par  = MX("p",3,1)
    qend,_,_,_ = integrator.call([q0,par])
    qe=MXFunction([q0,par],[qend])
    qe.init()
    qendJ=integrator.jacobian("x0","xf")
    qendJ.init()
    qendJ =qendJ.call([q0,par])[0]
    qeJ=MXFunction([q0,par],[qendJ])
    qeJ.init()

    qe.input(0).set(A)
    qe.input(1).set(B)
    qe.evaluate()
    print array(qe.output())

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

    t=ssym("t")
    q=ssym("y",2,1)
    p=ssym("p",1,1)
    # y
    # y'
    f=SXFunction(daeIn(x=q,p=p,t=t),daeOut(ode=vertcat([q[1],p[0]+q[1]**2 ])))
    f.init()
    
    integrator = CVodesIntegrator(f)
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
    q0   = MX("q0",2,1)
    par  = MX("p",1,1)
    qend,_,_,_ = integrator.call([q0,par])
    qe=MXFunction([q0,par],[qend])
    qe.init()
    qendJ=integrator.jacobian("x0","xf")
    qendJ.init()
    qendJ = qendJ.call([q0,par])[0]
    qeJ=MXFunction([q0,par],[qendJ])
    qeJ.init()

    qe.input(0).set(A)
    qe.input(1).set(p0)
    qe.evaluate()

    print qe.output()[0]
    print qe.output()[1]
    
    self.assertAlmostEqual(qe.output()[0],(2*y0-log(yc0**2/p0+1))/2-log(cos(arctan(yc0/sqrt(p0))+sqrt(p0)*te)),11,"Nonlin ODE")
    self.assertAlmostEqual(qe.output()[1],sqrt(p0)*tan(arctan(yc0/sqrt(p0))+sqrt(p0)*te),11,"Nonlin ODE")
    
    qeJ.input(0).set(A)
    qeJ.input(1).set(p0)
    qeJ.evaluate()
    
    Jr = array([[1,(sqrt(p0)*tan(sqrt(p0)*te+arctan(dy0/sqrt(p0)))-dy0)/(dy0**2+p0)],[0,(p0*tan(sqrt(p0)*te+arctan(dy0/sqrt(p0)))**2+p0)/(dy0**2+p0)]])
    self.checkarray(qeJ.output(),Jr,"jacobian of Nonlin ODE")
    
    #qe.setOption("ad_mode","reverse")
    #qe.init()
    Jf=qe.jacobian(0,0)
    Jf.init()
    Jf.input(0).set(A)
    Jf.input(1).set(p0)
    Jf.evaluate()
    self.checkarray(Jf.output(),Jr,"Jacobian of Nonlin ODE")
    
    #qe.setOption("ad_mode","forward")
    #qe.init();
    Jf=qe.jacobian(0,0)
    Jf.init()
    Jf.input(0).set(A)
    Jf.input(1).set(p0)
    Jf.evaluate()
    self.checkarray(Jf.output(),Jr,"Jacobian of Nonlin ODE")
    
    # Joel: This is no longer supported: might be a good idea to support again, though
    #qeJ=integrator.jac("x0","xf")
    #qeJ.init()
    #qeJ.input("x0").set(list(A)+[0,1,0,0])
    #qeJ.adjSeed("xf").set([0,0]+[0,1,0,0])
    #qeJ.evaluate(0,1)
    #print qeJ.output()
    #print qeJ.adjSens("x0")
    
    Jr = matrix([[(sqrt(p0)*(te*yc0**2-yc0+p0*te)*tan(arctan(yc0/sqrt(p0))+sqrt(p0)*te)+yc0**2)/(2*p0*yc0**2+2*p0**2)],[(sqrt(p0)*((te*yc0**2-yc0+p0*te)*tan(arctan(yc0/sqrt(p0))+sqrt(p0)*te)**2+te*yc0**2-yc0+p0*te)+(yc0**2+p0)*tan(arctan(yc0/sqrt(p0))+sqrt(p0)*te))/(sqrt(p0)*(2*yc0**2+2*p0))]])  
    
    #qe.setOption("ad_mode","reverse")
    #qe.init()
    Jf=qe.jacobian(1,0)
    Jf.init()
    Jf.input(0).set(A)
    Jf.input(1).set(p0)
    Jf.evaluate()
    self.checkarray(Jf.output(),Jr,"Jacobian of Nonlin ODE")
    
    #qe.setOption("ad_mode","forward")
    #qe.init()
    Jf=qe.jacobian(1,0)
    Jf.init()
    Jf.input(0).set(A)
    Jf.input(1).set(p0)
    Jf.evaluate()
    self.checkarray(Jf.output(),Jr,"Jacobian of Nonlin ODE")
    
    qendJ=integrator.jacobian("p","xf")
    qendJ.init()
    qendJ = qendJ.call([q0,par])[0]
    qeJ=MXFunction([q0,par],[qendJ])
    qeJ.init()

    qeJ.input(0).set(A)
    qeJ.input(1).set(p0)
    qeJ.evaluate()
    
    self.checkarray(qeJ.output(),Jr,"jacobian of Nonlin ODE")
    
    
    
    
    qeJf=MXFunction([q0,par],[vec(qeJ.call([q0,par])[0])])
    #qeJf.setOption("ad_mode","reverse")
    qeJf.init()
    
    H=qeJf.jacobian(0,0)
    H.init()
    H.input(0).set(A)
    H.input(1).set(p0)
    H.evaluate()
    def sec(x):
      return 1.0/cos(x)
    Hr = array([[0,0],[0,-(2*yc0*tan(arctan(yc0)+te))/(yc0**4+2*yc0**2+1)+sec(arctan(yc0)+te)**2/(yc0**4+2*yc0**2+1)+(2*yc0**2)/(yc0**4+2*yc0**2+1)-1/(yc0**2+1)],[0,0],[0,-(2*yc0*tan(arctan(yc0)+te)**2)/(yc0**4+2*yc0**2+1)+(2*sec(arctan(yc0)+te)**2*tan(arctan(yc0)+te))/(yc0**4+2*yc0**2+1)-(2*yc0)/(yc0**4+2*yc0**2+1)]])
    print array(H.output())
    print Hr
    
    # Joel: As above, this is no longer supported
    #qeJ=integrator.jac("x0","xf")
    #qeJ.init()
    #qeJ.input("x0").set(list(A)+[0,1,0,0])
    #qeJ.adjSeed("xf").set([0,0]+[0,1,0,0])
    #qeJ.evaluate(0,1)
    #print qeJ.output()
    #print qeJ.adjSens("x0")
    

  def test_hessian2D(self):
    self.message("hessian")
    N=2

    x0_ = DMatrix([1,0.1])
    A_  = DMatrix([[3,1],[0.74,4]])

    A = ssym("A",N,N)
    x = ssym("x",N)

    ode = SXFunction(daeIn(x=x, p=vec(A)),daeOut(ode=mul(A,x)))
    I = CVodesIntegrator(ode)
    I.setOption("fsens_err_con", True)
    I.setOption('reltol',1e-12)
    I.init()
    I.input("x0").set(x0_)
    I.input("p").set(vec(A_))
    I.evaluate()

    q0=MX("q0",N)
    p=MX("p",N*N)
    qe = MXFunction([q0,p],I.call([q0,p]))
    qe.init()

    JT = MXFunction([q0,p],[qe.jac(1,0).T])
    JT.init()

    H  = JT.jacobian(1)
    H.init()
    H.input(0).set(x0_)
    H.input(1).set(vec(A_))
    H.evaluate()

    H1 = DMatrix(H.output())
    
    ## Joel: Only Hessians of scalar functions allowed
    #H = qe.hessian(1)
    #H.init()
    #H.input(0).set(x0_)
    #H.input(1).set(vec(A_))
    #H.evaluate()
    #H2 = DMatrix(H.output())
    
    #self.checkarray(H1,H2,"hessian")
    
  def test_issue535(self):
    self.message("regression test for #535")
    t=ssym("t")
    x=ssym("x")
    rx=ssym("rx")
    p=ssym("p")
    dp=ssym("dp")

    z=ssym("z")
    rz=ssym("rz")
    rp=ssym("rp")
    f = SXFunction(daeIn(**{'x': x, 'z': z}),daeOut(**{'alg': x-z, 'ode': z}))
    f.init()
    g = SXFunction(rdaeIn(**{'x': x, 'z': z, 'rx': rx, 'rz': rz}),rdaeOut(**{'alg': x-rz, 'ode': rz}))
    g.init()

    integrator = IdasIntegrator(f,g)
    integrator.setOption({'calc_ic': True, 'tf': 2.3, 'reltol': 1e-10, 'augmented_options': {'reltol': 1e-09, 'abstol': 1e-09 }, 'calc_icB': True, 'abstol': 1e-10, 't0': 0.2})
    integrator.init()

    integrator.input("x0").set(7.1)
    if not integrator.input("p").empty():
      integrator.input("p").set(2)
    if not integrator.input("rx0").empty():
      integrator.input("rx0").set(0.13)
    if not integrator.input("rp").empty():
      integrator.input("rp").set(0.127)

    integrator.evaluate(0,0)

    integrator.fwdSeed(0).set([1])
    integrator.evaluate(1,0) # fail
    
if __name__ == '__main__':
    unittest.main()

