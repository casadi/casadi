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
import numpy
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
  load_integrator("cvodes")
  integrators.append(("cvodes",["ode"],{"abstol": 1e-15,"reltol":1e-15,"fsens_err_con": True,"quad_err_con": False}))
except:
  pass

try:
  load_integrator("idas")
  integrators.append(("idas",["dae","ode"],{"abstol": 1e-15,"reltol":1e-15,"fsens_err_con": True,"calc_icB":True}))
except:
  pass

integrators.append(("collocation",["dae","ode"],{"rootfinder":"kinsol","number_of_finite_elements": 18}))

integrators.append(("rk",["ode"],{"number_of_finite_elements": 1000}))

print("Will test these integrators:")
for cl, t, options in integrators:
  print(cl, " : ", t)

class Integrationtests(casadiTestCase):

  @slow()
  def test_full(self):
    num = self.num
    tc = DM(n.linspace(0,num['tend'],100))

    t=SX.sym("t")
    q=SX.sym("q")
    p=SX.sym("p")

    dae = {'t':t, 'x':q, 'p':p, 'ode':q/p*t**2}
    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    opts["fsens_err_con"] = True
    #opts["verbose"] = True
    opts["t0"] = 0
    opts["tf"] = 2.3
    integrator = casadi.integrator("integrator", "cvodes", dae, opts)
    tf = 2.3

    solution = Function("solution", {'x0':q, 'p':p, 'xf':q*exp(tf**3/(3*p))},
                        casadi.integrator_in(), casadi.integrator_out())
    f_in = {}
    f_in["x0"]=0.3
    f_in["p"]=0.7

    self.checkfunction(integrator,solution,inputs=f_in,digits=6)

  def test_tools_trivial(self):
    num = self.num

    x = SX.sym("x")
    p = SX.sym("p",0)
    tf = SX.sym("tf")
    f=Function("f", [x,p], [x])

    for integrator in [
         simpleRK(f),
         simpleIRK(f),
         simpleIntegrator(f)
       ]:

      solution = Function("solution", [x,p,tf], [x*exp(tf)], ["x0","p","h"], ["xf"])

      f_in = {"x0":1,"h":1}
      self.checkfunction(integrator,solution,inputs=f_in,digits=3)

  @slow()
  def test_tools(self):
    num = self.num

    t=SX.sym("t")
    q=SX.sym("q")
    p=SX.sym("p")
    t0=SX.sym("t0")
    q0=SX.sym("q0")
    f=Function("f", [vertcat(*(q,t)),p],[vertcat(*(q/p*t**2,1))])
    for integrator in [
            simpleRK(f,500),
            simpleIRK(f,50),
            simpleIntegrator(f)
            ]:

      solution = Function('solver', {"x0": vertcat(*(q0,t0)),"p":p,"h":t, "xf":vertcat(*[q0*exp(((t0+t)**3-t0**3)/(3*p)),t0+t])}, integrator.name_in(), integrator.name_out())

      f_in = {}
      f_in["x0"]=DM([0.3,0])
      f_in["p"]=0.7
      f_in["h"]=1

      self.checkfunction(integrator,solution,inputs=f_in,digits=3)

  @memory_heavy()
  def test_jac(self):
    self.message("Test exact jacobian #536")
    # This test is not automized, but works by inspection only.
    # To activate, recompile after ucnommenting the printout lines in cvodes.c, near "Used for validating casadi#536"
    #return
    DM.set_precision(18)

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
          x = SX.sym("x")
          dummyIntegrator = casadi.integrator("dummyIntegrator", Integrator, {'x':x, 'ode':x})
          if p_features[0] in features:

            dae = din.copy()
            dae.update(dout)
            dae.update(rdin)
            for k, v in rdout.items(): dae['r'+k] = v

            for k in list(solution.keys()):
              solution[k] = substitute(solution[k],vertcat(*[tstart,tend]),vertcat(*[tstart_,tend_]))

            sol = solutionin.copy()
            sol.update(solution)
            fs = Function("fs", sol, casadi.integrator_in(), casadi.integrator_out())

            def itoptions():
                yield {"newton_scheme": "direct"}
                yield {"newton_scheme": "gmres"}
                yield {"newton_scheme": "bcgstab"}
                # yield {"newton_scheme": "tfqmr"} # Bug in Sundials?

            for f_options in itoptions():
                message = "f_options: %s" % str(f_options)
                print(message)
                opts = {}
                opts["t0"] = tstart_
                opts["tf"] = tend_
                for op in (options, f_options):
                  for (k,v) in list(op.items()):
                    opts[k] = v
                integrator = casadi.integrator("integrator", Integrator, dae, opts)
                integrator_in = {}
                for k,v in list(point.items()):
                  if not fs.sparsity_in(k).is_empty():
                    integrator_in[k]=v

                integrator_out = integrator(**integrator_in)
                fs_out = fs(**integrator_in)
                print("res=",integrator_out["xf"]-fs_out["xf"], fs_out["xf"])
                print("Rres=",integrator_out["rxf"]-fs_out["rxf"], fs_out["rxf"])
                # self.checkarray(integrator_out["rxf"],fs_out["rxf"],digits=4)
                stats = integrator.stats()

                print(stats)
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
      t=SX.sym("t")
      x=SX.sym("x")
      rx=SX.sym("rx")
      p=SX.sym("p")
      dp=SX.sym("dp")

      z=SX.sym("z")
      rz=SX.sym("rz")
      rp=SX.sym("rp")
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
      yield (["ode"],{'x':q,'p':p},{'ode': vertcat(*[q[1],p[0]+q[1]**2 ])},{},{},{'x0':q, 'p': p} ,{'xf': vertcat(*[s1,s2])},{'x0': A, 'p': p0},(0,0.4) )

    for p_features, din, dout, rdin, rdout, solutionin, solution, point, (tstart_, tend_) in checks():

      for Integrator, features, options in integrators:
        self.message(Integrator)
        x = SX.sym("x")
        dummyIntegrator = casadi.integrator("dummyIntegrator", Integrator, {'x':x, 'ode':x})
        if p_features[0] in features:

          dae = din.copy()
          dae.update(dout)
          dae.update(rdin)
          for k, v in rdout.items(): dae['r'+k] = v

          for k in list(solution.keys()):
            solution[k] = substitute(solution[k],vertcat(*[tstart,tend]),vertcat(*[tstart_,tend_]))

          sol = solutionin.copy()
          sol.update(solution)
          fs = Function("fs", sol, casadi.integrator_in(), casadi.integrator_out())

          def itoptions():
            yield {"newton_scheme": "direct"}
            yield {"newton_scheme": "gmres"}
            yield {"newton_scheme": "bcgstab"}
            #yield {"newton_scheme": "tfqmr"} # Bug in Sundials?

          for f_options in itoptions():
              message = "f_options: %s" % str(f_options)
              print(message)

              opts = {}
              opts["t0"] = tstart_
              opts["tf"] = tend_
              for op in (options, f_options):
                 for (k,v) in list(op.items()):
                    opts[k] = v
              integrator = casadi.integrator("integrator", Integrator, dae, opts)
              integrator_in = {}

              for k,v in list(point.items()):
                if not integrator.sparsity_in(k).is_empty():
                  integrator_in[k]=v

              self.checkfunction(integrator,fs,inputs=integrator_in,gradient=False,hessian=False,sens_der=False,evals=False,digits=4,digits_sens=4,failmessage=message,verbose=False)




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
        yield (["ode"],{'x':x},{'ode':x},{'x':x,'rx':rx},{'ode':SX(0)},si,{'rxf': rx},pointA,ti)
        yield (["ode"],{'x':x},{'ode':x},{'x':x,'rx':rx},{'ode':SX(1)},si,{'rxf': rx+tend-tstart},pointA,ti)
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
        yield (["ode"],{'x':q,'p':p},{'ode': vertcat(*[q[1],p[0]+q[1]**2 ])},{},{},{'x0':q, 'p': p} ,{'xf': vertcat(*[s1,s2])},{'x0': A, 'p': p0},(0,0.4) )

      for tt in checks():
        print(tt)
        for p_features, din, dout, rdin, rdout, solutionin, solution, point, (tstart_, tend_) in variations(*tt):
          if p_features[0] in features:
            message = "%s: %s => %s, %s => %s, explicit (%s) tstart = %f" % (Integrator,str(din),str(dout),str(rdin),str(rdout),str(solution),tstart_)
            print(message)

            dae = din.copy()
            dae.update(dout)
            dae.update(rdin)
            for k, v in rdout.items(): dae['r'+k] = v

            for k in list(solution.keys()):
              solution[k] = substitute(solution[k],vertcat(*[tstart,tend]),vertcat(*[tstart_,tend_]))

            sol = solutionin.copy()
            sol.update(solution)
            fs = Function("fs", sol, casadi.integrator_in(), casadi.integrator_out())

            opts = dict(options)
            opts["t0"] = tstart_
            if Integrator in ('cvodes', 'idas'):
              opts["abstol"] = 1e-9
              opts["reltol"] = 1e-9
            opts["tf"] = tend_
            if Integrator=='idas':
              opts["init_xdot"] = DM(point["x0"]).nonzeros()
              opts["calc_icB"] = True
              opts["augmented_options"] = {"init_xdot":None, "abstol":1e-9,"reltol":1e-9}
            integrator = casadi.integrator("integrator", Integrator, dae, opts)
            integrator_in = {}

            for k,v in list(point.items()):
              if not integrator.sparsity_in(k).is_empty():
                integrator_in[k]=v

            self.checkfunction(integrator,fs,inputs=integrator_in,gradient=False,hessian=False,sens_der=False,evals=False,digits=4,digits_sens=4,failmessage=message,verbose=False)


  def setUp(self):
    # Reference solution is x0 e^((t^3-t0^3)/(3 p))
    t=SX.sym("t")
    x=SX.sym("x")
    p=SX.sym("p")
    dae={'t':t, 'x':x, 'p':p, 'ode':x/p*t**2}
    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    opts["fsens_err_con"] = True
    #opts["verbose"] = True
    opts["t0"] = 0
    opts["tf"] = 2.3
    integrator = casadi.integrator("integrator", "cvodes", dae, opts)
    q0   = MX.sym("q0")
    par  = MX.sym("p")

    qend = integrator(x0=q0, p=par)["xf"]

    qe=Function("qe", [q0,par],[qend])
    self.integrator = integrator
    self.qe=qe
    self.qend=qend
    self.q0=q0
    self.par=par
    self.dae = dae
    self.num={'tend':2.3,'q0':7.1,'p':2}
    pass

  def test_eval2(self):
    self.message('CVodes integration: evaluation with Function indirection')
    num=self.num
    qend=self.qend

    par=self.par
    q0=self.q0
    qe=Function("qe", [q0,par],[qend[0]])

    f = Function("f", [q0], [qe(q0,MX(num['p']))])
    f_out = f([num['q0']])
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(f_out[0][0],q0*exp(tend**3/(3*p)),9,"Evaluation output mismatch")

  def test_issue92b(self):
    self.message("regression check for issue 92")
    t=SX.sym("t")
    x=SX.sym("x")
    y=SX.sym("y")
    dae = {'t':t, 'x':vertcat(*[x,y]), 'ode':vertcat(*[x,(1+1e-9)*x])}
    opts = {}
    opts["fsens_err_con"] = True
    opts["t0"] = 0
    opts["tf"] = 1
    integrator = casadi.integrator("integrator", "cvodes", dae, opts)
    integrator_in = {}

    # Pass inputs
    integrator_in["x0"]=DM([1,0])
    ## Integrate
    integrator_out = integrator(**integrator_in)
    # print result
    print(integrator_out["xf"])

  def test_issue92(self):
    self.message("regression check for issue 92")
    t=SX.sym("t")
    x=SX.sym("x")
    var = MX.sym("var",2,1)

    q = vertcat(*[x,SX.sym("problem")])

    dq=vertcat(*[x,x])
    dae = {'t':t, 'x':q, 'ode':dq}
    opts = {}
    opts["fsens_err_con"] = True
    opts["reltol"] = 1e-12
    opts["t0"] = 0
    opts["tf"] = 1
    integrator = casadi.integrator("integrator", "cvodes", dae, opts)

    qend = integrator(x0=var)["xf"]

    f = Function("f", [var],[qend[0]])

    J=f.jacobian_old(0, 0)
    J_out = J([1,0])
    print("jac=",J_out[0].nz[0]-exp(1))
    self.assertAlmostEqual(J_out[0][0,0],exp(1),5,"Evaluation output mismatch")

  def test_eval(self):
    self.message('CVodes integration: evaluation')
    num=self.num
    qe=self.qe
    qe_out = qe([num['q0']], [num['p']])
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(qe_out[0][0],q0*exp(tend**3/(3*p)),9,"Evaluation output mismatch")


  def test_jac1(self):
    self.message('CVodes integration: jacobian to q0')
    num=self.num
    J=self.qe.jacobian_old(0, 0)
    J_out = J([num['q0']], [num['p']])
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(J_out[0][0],exp(tend**3/(3*p)),9,"Evaluation output mismatch")

  def test_jac2(self):
    self.message('CVodes integration: jacobian to p')
    num=self.num
    J=self.qe.jacobian_old(1, 0)
    J_out = J([num['q0']], [num['p']])
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(J_out[0][0],-(q0*tend**3*exp(tend**3/(3*p)))/(3*p**2),9,"Evaluation output mismatch")

  def test_bug_repeat(self):
    num={'tend':2.3,'q0':[0,7.1,7.1],'p':2}
    self.message("Bug that appears when rhs contains repeats")
    A=array([1,0.1,1])
    p0 = 1.13
    y0=A[0]
    yc0=dy0=A[1]
    te=0.4

    t=SX.sym("t")
    q=SX.sym("y",3,1)
    p=SX.sym("p")

    dh = p+q[0]**2
    dae = {'x':q, 'p':p, 't':t, 'ode':vertcat(*[dh ,q[0],dh])}

    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    #opts["verbose"] = True
    opts["fsens_err_con"] = True
    opts["steps_per_checkpoint"] = 10000
    opts["t0"] = 0
    opts["tf"] = te
    integrator = casadi.integrator("integrator", "cvodes", dae, opts)

    q0   = MX.sym("q0",3,1)
    par  = MX.sym("p",1,1)
    qend = integrator(x0=q0, p=par)["xf"]
    qe=Function("qe", [q0,par],[qend])

    #J=self.qe.jacobian_old(2, 0)
    J=qe.jacobian_old(0, 0)
    J_out = J(A, p0)
    outA=J_out[0].full()
    dae={'x':q, 'p':p, 't':t, 'ode':vertcat(*[dh ,q[0],(1+1e-9)*dh])}

    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    #opts["verbose"] = True
    opts["fsens_err_con"] = True
    opts["steps_per_checkpoint"] = 10000
    opts["t0"] = 0
    opts["tf"] = te
    integrator = casadi.integrator("integrator", "cvodes", dae, opts)

    q0   = MX.sym("q0",3,1)
    par  = MX.sym("p",1,1)
    qend = integrator(x0=q0, p=par)["xf"]
    qe=Function("qe", [q0,par],[qend])

    #J=self.qe.jacobian_old(2)
    J=qe.jacobian_old(0, 0)
    J_out = J(A, p0)
    outB=J_out[0].full()
    print(outA-outB)

  def test_hess3(self):
    self.message('CVodes integration: hessian to p: Jacobian of integrator.jacobian')
    num=self.num
    J=self.integrator.jacobian_old(self.integrator.index_in("p"),self.integrator.index_out("xf"))
    H=J.jacobian_old(J.index_in("p"), 0)
    H_in = {}
    H_in["x0"]=num['q0']
    H_in["p"]=num['p']
    H_out = H(**H_in)
    num=self.num
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(H_out["jac_jac_xf_p_p"][0],(q0*tend**6*exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*exp(tend**3/(3*p)))/(3*p**3),9,"Evaluation output mismatch")

  def test_hess4(self):
    self.message('CVodes integration: hessian to p: Jacobian of integrator.jacobian indirect')
    num=self.num
    J=self.integrator.jacobian_old(self.integrator.index_in("p"),self.integrator.index_out("xf"))

    q0=MX.sym("q0")
    p=MX.sym("p")
    Ji = J(x0=q0,p=p)
    Ji["q0"] = q0
    Ji["p"] = p
    Ji = Function("Ji", Ji, ["q0", "p"], J.name_out())
    H=Ji.jacobian_old(1, 0)
    H_out = H([num['q0']], [num['p']])
    num=self.num
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(H_out[0][0],(q0*tend**6*exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*exp(tend**3/(3*p)))/(3*p**3),9,"Evaluation output mismatch")

  def test_hess5(self):
    self.message('CVodes integration: hessian to p in an MX tree')
    num=self.num
    q0=MX.sym("q0")
    p=MX.sym("p")

    sol = self.integrator(x0=q0,p=p)
    sol["q0"] = q0
    sol["p"] = p
    JT = Function("JT", [q0,p],[jacobian(sol['xf'], sol['p']).T])
    JT_out = JT([num['q0']], [num['p']])
    print(JT_out)

    H  = JT.jacobian_old(1, 0)
    H_out = H([num['q0']], [num['p']])
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(H_out[0],(q0*tend**6*exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*exp(tend**3/(3*p)))/(3*p**3),9,"Evaluation output mismatch")

  def test_hess6(self):
    self.message('CVodes integration: hessian to p in an MX tree')
    num=self.num
    q0=MX.sym("q0")
    p=MX.sym("p")

    sol = self.integrator(x0=q0,p=p)
    sol["q0"] = q0
    sol["p"] = p
    qe = Function("qe", sol, ["q0", "p"], casadi.integrator_out())

    H = qe.hessian_old(1, 0)
    H_out = H([num['q0']], [num['p']])
    num=self.num
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(H_out[0][0],(q0*tend**6*exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*exp(tend**3/(3*p)))/(3*p**3),8,"Evaluation output mismatch")

  def test_glibcbug(self):
    self.message("former glibc error")
    A=array([2.3,4.3,7.6])
    B=array([[1,2.3,4],[-2,1.3,4.7],[-2,6,9]])

    te=0.7
    t=SX.sym("t")
    q=SX.sym("q",3,1)
    p=SX.sym("p",9,1)
    dae = {'t':t, 'x':q, 'p':p, 'ode':mtimes(c.reshape(p,3,3),q)}
    opts = {}
    opts["fsens_err_con"] = True
    opts["steps_per_checkpoint"] = 1000
    opts["t0"] = 0
    opts["tf"] = te
    integrator = casadi.integrator("integrator", "cvodes", dae, opts)
    q0   = MX.sym("q0",3,1)
    par  = MX.sym("p",9,1)
    qend = integrator(x0=q0, p=par)["xf"]
    qe=integrator.jacobian_old(integrator.index_in("p"),integrator.index_out("xf"))
    qe = qe(x0=q0,p=par)['jac_xf_p']
    qef=Function("qef", [q0,par],[qe])
    qef_out = qef(A, B.ravel())

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

    dae = {'t':t, 'x':q, 'p':p, 'ode':mtimes(c.reshape(p,3,3),q)}
    opts = {}
    opts["fsens_err_con"] = True
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    #opts["verbose"] = True
    opts["steps_per_checkpoint"] = 10000
    opts["t0"] = 0
    opts["tf"] = te
    integrator = casadi.integrator("integrator", "cvodes", dae, opts)

    q0   = MX.sym("q0",3,1)
    par  = MX.sym("p",9,1)
    qend = integrator(x0=q0, p=par)["xf"]
    qe=Function("qe", [q0,par],[qend])
    qendJ=integrator.jacobian_old(integrator.index_in("x0"),integrator.index_out("xf"))
    qendJ = qendJ(x0=q0,p=par)['jac_xf_x0']

    qeJ=Function("qeJ", [q0,par],[qendJ])

    qendJ2=integrator.jacobian_old(integrator.index_in("x0"),integrator.index_out("xf"))
    qendJ2 = qendJ2(x0=q0,p=par)['jac_xf_x0']

    qeJ2=Function("qeJ2", [q0,par],[qendJ2])
    qe_out = qe(A, vec(B))
    self.checkarray(numpy.dot(Be,A)/1e3,qe_out/1e3,"jacobian('x0','xf')")
    qeJ_out = qeJ(A, vec(B))
    self.checkarray(qeJ_out/1e3,Be/1e3,"jacobian('x0','xf')")

    qeJ2_out = qeJ2(A, vec(B))

    return # this should return identical zero
    H=qeJ.jacobian_old(0, 0)
    H_out = H(A, vec(B))
    print(array(H_out[0]))


  def test_mathieu_system(self):
    self.message("Mathieu ODE")
    A=array([0.3,1.2])
    B=array([1.3,4.3,2.7])
    te=0.7

    t=SX.sym("t")
    q=SX.sym("q",2,1)
    p=SX.sym("p",3,1)

    dae = {'x':q, 'p':p, 't':t, 'ode':vertcat(*[q[1],(p[0]-2*p[1]*cos(2*p[2]))*q[0]])}
    opts = {}
    opts["fsens_err_con"] = True
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    #opts["verbose"] = True
    opts["steps_per_checkpoint"] = 10000
    opts["t0"] = 0
    opts["tf"] = te
    integrator = casadi.integrator("integrator", "cvodes", dae, opts)

    q0   = MX.sym("q0",2,1)
    par  = MX.sym("p",3,1)
    qend = integrator(x0=q0, p=par)["xf"]
    qe=Function("qe", [q0,par],[qend])
    qendJ=integrator.jacobian_old(integrator.index_in("x0"), integrator.index_out("xf"))
    qendJ =qendJ(x0=q0,p=par)['jac_xf_x0']
    qeJ=Function("qeJ", [q0,par],[qendJ])
    qe_out = qe(A, B)
    print(array(qe_out))

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
    dae = {'x':q, 'p':p, 't':t, 'ode':vertcat(*[q[1],p[0]+q[1]**2])}
    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    #opts["verbose"] = True
    opts["steps_per_checkpoint"] = 10000
    opts["fsens_err_con"] = True
    opts["t0"] = 0
    opts["tf"] = te
    integrator = casadi.integrator("integrator", "cvodes", dae, opts)

    t0   = MX(0)
    tend = MX(te)
    q0   = MX.sym("q0",2,1)
    par  = MX.sym("p",1,1)
    qend = integrator(x0=q0, p=par)["xf"]
    qe=Function("qe", [q0,par],[qend])
    qendJ=integrator.jacobian_old(integrator.index_in("x0"), integrator.index_out("xf"))
    qendJ = qendJ(x0=q0, p=par)['jac_xf_x0']
    qeJ=Function("qeJ", [q0,par],[qendJ])
    qe_out = qe(A, p0)

    print(qe_out[0])
    print(qe_out[1])

    self.assertAlmostEqual(qe_out[0],(2*y0-log(yc0**2/p0+1))/2-log(cos(arctan(yc0/sqrt(p0))+sqrt(p0)*te)),11,"Nonlin ODE")
    self.assertAlmostEqual(qe_out[1],sqrt(p0)*tan(arctan(yc0/sqrt(p0))+sqrt(p0)*te),11,"Nonlin ODE")

    qeJ_out = qeJ(A, p0)

    Jr = array([[1,(sqrt(p0)*tan(sqrt(p0)*te+arctan(dy0/sqrt(p0)))-dy0)/(dy0**2+p0)],[0,(p0*tan(sqrt(p0)*te+arctan(dy0/sqrt(p0)))**2+p0)/(dy0**2+p0)]])
    self.checkarray(qeJ_out,Jr,"jacobian of Nonlin ODE")

    Jf=qe.jacobian_old(0,0)
    Jf_out = Jf(A, p0)
    self.checkarray(Jf_out[0],Jr,"Jacobian of Nonlin ODE")

    Jf=qe.jacobian_old(0,0)
    Jf_out = Jf(A, p0)
    self.checkarray(Jf_out[0],Jr,"Jacobian of Nonlin ODE")

    Jr = numpy.matrix([[(sqrt(p0)*(te*yc0**2-yc0+p0*te)*tan(arctan(yc0/sqrt(p0))+sqrt(p0)*te)+yc0**2)/(2*p0*yc0**2+2*p0**2)],[(sqrt(p0)*((te*yc0**2-yc0+p0*te)*tan(arctan(yc0/sqrt(p0))+sqrt(p0)*te)**2+te*yc0**2-yc0+p0*te)+(yc0**2+p0)*tan(arctan(yc0/sqrt(p0))+sqrt(p0)*te))/(sqrt(p0)*(2*yc0**2+2*p0))]])

    Jf=qe.jacobian_old(1,0)
    Jf_out = Jf(A, p0)
    self.checkarray(Jf_out[0],Jr,"Jacobian of Nonlin ODE")
    Jf=qe.jacobian_old(1,0)
    Jf_out = Jf(A, p0)
    self.checkarray(Jf_out[0],Jr,"Jacobian of Nonlin ODE")

    qendJ=integrator.jacobian_old(integrator.index_in("p"),integrator.index_out("xf"))
    qendJ = qendJ(x0=q0,p=par)['jac_xf_p']
    qeJ=Function("qeJ", [q0,par],[qendJ])

    qeJ_out = qeJ(A, p0)

    self.checkarray(qeJ_out,Jr,"jacobian of Nonlin ODE")

    qeJf=Function("qeJf", [q0,par],[vec(qeJ(q0,par))])

    H=qeJf.jacobian_old(0,0)
    H_out = H(A, p0)
    def sec(x):
      return 1.0/cos(x)
    Hr = array([[0,0],[0,-(2*yc0*tan(arctan(yc0)+te))/(yc0**4+2*yc0**2+1)+sec(arctan(yc0)+te)**2/(yc0**4+2*yc0**2+1)+(2*yc0**2)/(yc0**4+2*yc0**2+1)-1/(yc0**2+1)],[0,0],[0,-(2*yc0*tan(arctan(yc0)+te)**2)/(yc0**4+2*yc0**2+1)+(2*sec(arctan(yc0)+te)**2*tan(arctan(yc0)+te))/(yc0**4+2*yc0**2+1)-(2*yc0)/(yc0**4+2*yc0**2+1)]])
    print(array(H_out[0]))
    print(Hr)

  def test_missing_symbols(self):
    x = SX.sym('x'); z = SX.sym('z'); p = SX.sym('p')
    dae = {'x':x, 'z':z, 'ode':z+p, 'alg':z*cos(z)-x} # p forgotten here
    with self.assertInException("[p] are free"):
      integrator('F', 'idas', dae)
  def test_hessian2D(self):
    self.message("hessian")
    N=2

    x0_ = DM([1,0.1])
    A_  = DM([[3,1],[0.74,4]])

    A = SX.sym("A",N,N)
    x = SX.sym("x",N)

    dae = {'x':x, 'p':vec(A), 'ode':mtimes(A,x)}
    I = casadi.integrator("I", "cvodes", dae, {"fsens_err_con": True, 'reltol' : 1e-12})
    I_in =  {}
    I_in["x0"]=x0_
    I_in["p"]=vec(A_)
    I_out = I(**I_in)

    q0=MX.sym("q0",N)
    p=MX.sym("p",N*N)

    sol = I(x0=q0,p=p)
    sol["q0"] = q0
    sol["p"] = p
    qe = Function("qe", sol, ["q0", "p"], casadi.integrator_out())

    JT = Function("JT", [q0,p],[jacobian(sol['xf'], sol['p']).T])

    H  = JT.jacobian_old(1, 0)
    H_out = H(x0_, vec(A_))

    H1 = H_out[0]

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
    dae = {'x': x, 'z': z, 'rx': rx, 'rz': rz, 'alg': x-z, 'ode': z, 'ralg': x-rz, 'rode': rz}

    integrator = casadi.integrator("integrator", "idas", dae, {'calc_ic': True, 'tf': 2.3, 'reltol': 1e-10, 'augmented_options': {'reltol': 1e-09, 'abstol': 1e-09 }, 'calc_icB': True, 'abstol': 1e-10, 't0': 0.2})
    integrator_in = {}

    integrator_in["x0"]=7.1
    if not integrator.sparsity_in("p").is_empty():
      integrator_in["p"]=2
    if not integrator.sparsity_in("rx0").is_empty():
      integrator_in["rx0"]=0.13
    if not integrator.sparsity_in("rp").is_empty():
      integrator_in["rp"]=0.127

    integrator_out = integrator(**integrator_in)

  def test_collocationPoints(self):
    self.message("collocation points")
    with self.assertRaises(Exception):
      collocation_points(0,"radau")
    with self.assertRaises(Exception):
      collocation_points(10,"radau")
    with self.assertRaises(Exception):
      collocation_points(0,"legendre")
    with self.assertRaises(Exception):
      collocation_points(10,"legendre")
    with self.assertRaises(Exception):
      collocation_points(1,"foo")

    for k in range(1,10):
      r = [0] + collocation_points(k,"radau")
      self.assertEqual(len(r),k+1)
      self.checkarray(DM(r[-1]),DM([1]))
    for k in range(1,10):
      r = [0] + collocation_points(k,"legendre")
      self.assertEqual(len(r),k+1)


  @memory_heavy()
  def test_thread_safety(self):
    x = MX.sym('x')
    rx = MX.sym('rx')
    t = MX.sym('t')
    p = MX.sym('p')
    for Integrator, features, options in integrators:
      intg = integrator("Integrator",Integrator,{"x":x,"ode":-x}, options)

      intg_par = intg.map(40, 'thread',4)
      res = intg_par(x0=numpy.linspace(0, 10, 40))
      self.checkarray(norm_inf(res["xf"].T-exp(-1)*numpy.linspace(0, 10, 40)),0, digits=5)

      if Integrator=="cvodes": continue
      if Integrator=="idas": continue

      intg = integrator("Integrator",Integrator,{"x":x,"rx":rx,"ode":-x,"rode": rx}, options)

      intg_par = intg.map(40, 'thread',2)
      res = intg_par(rx0=numpy.linspace(0, 10, 40), x0=numpy.linspace(0, 10, 40))
      self.checkarray(norm_inf(res["xf"].T-exp(-1)*numpy.linspace(0, 10, 40)),0, digits=5)
      self.checkarray(norm_inf(res["rxf"].T-exp(1)*numpy.linspace(0, 10, 40)),0, digits=5)

if __name__ == '__main__':
    unittest.main()
