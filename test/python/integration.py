#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             KU Leuven. All rights reserved.
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

import casadi as ca
import numpy as np
import casadi as c
import numpy
import numpy as n
import unittest
from types import *
from numpy import array
from helpers import args, casadiTestCase, hessian_old, jacobian_old, memory_heavy, requires_integrator, requires_nlpsol, slow
import copy

scipy_available = True
try:
  import scipy.special
  from scipy.linalg import expm
except:
  scipy_available = False

integrators = []

try:
  ca.load_integrator("cvodes")
  integrators.append(("cvodes",["ode"],{"abstol": 1e-15,"reltol":1e-15,"fsens_err_con": True,"quad_err_con": False}))
except:
  pass

try:
  ca.load_integrator("idas")
  integrators.append(("idas",["dae","ode"],{"abstol": 1e-15,"reltol":1e-15,"fsens_err_con": True,"calc_icB":True}))
except:
  pass

integrators.append(("collocation",["dae","ode"],{"rootfinder":"newton","number_of_finite_elements": 18}))

integrators.append(("rk",["ode"],{"number_of_finite_elements": 1000}))

integrators.append(("collocation",["dae","ode"],{"rootfinder":"newton","number_of_finite_elements": 18,"simplify":True,"rootfinder":"fast_newton"}))

integrators.append(("rk",["ode"],{"number_of_finite_elements": 1000,"simplify":True}))

integrators.append(("rosenbrock",["dae","ode"],{"abstol": 1e-12,"reltol":1e-12}))

# Integrators with zero-crossing events (option transition), by feature: events at all, of
# DAEs (also sensitivities), adjoints, code generation; digits of the results. The event
# handling of Integrator (cvodes, idas, rk, collocation) gets forward sensitivities through
# events of a DAE wrong: the rates after the transition use inconsistent algebraic states
event_integrators = []

for plugin in ["cvodes", "idas"]:
  try:
    ca.load_integrator(plugin)
  except:
    continue
  event_integrators.append((plugin, {"abstol": 1e-12, "reltol": 1e-12, "quad_err_con": True,
    "max_event_iter": 20}, {"events": True, "events_dae": False,
    "events_adjoint": False, "codegen": False, "digits": 5}))

event_integrators.append(("rk", {"number_of_finite_elements": 100, "max_event_iter": 20},
  {"events": True, "events_dae": False, "events_adjoint": False, "codegen": False, "digits": 7}))

event_integrators.append(("collocation", {"number_of_finite_elements": 20, "max_event_iter": 20},
  {"events": True, "events_dae": False, "events_adjoint": False, "codegen": False, "digits": 7}))

for plugin, options in [("erk", {}), ("erk", {"method": "rk4", "step_size": 0.002}),
                        ("rosenbrock", {}), ("rosenbrock", {"step_size": 0.002})]:
  event_integrators.append((plugin, dict(options, abstol=1e-11, reltol=1e-11, max_tape_steps=600),
    {"events": True, "events_dae": plugin == "rosenbrock", "events_adjoint": True,
     "codegen": True, "digits": 8}))

def bounce_transition(nz=0):
  # Bouncing ball: the velocity is reversed and scaled by the restitution p; z passes through
  x, z, p = ca.SX.sym("x", 2), ca.SX.sym("z", nz), ca.SX.sym("p")
  index, t, u = ca.SX.sym("index"), ca.SX.sym("t"), ca.SX(0, 1)
  return ca.Function("transition", [index, t, x, z, p, u], [ca.vertcat(x[0], -p*x[1]), z],
    ["index", "t", "x", "z", "p", "u"], ["post_x", "post_z"])

def fd_jacobian(G, v0, h=1e-6):
  # Central differences of the first output of G at v0
  v0 = ca.DM(v0)
  f = lambda v: (lambda r: r[0] if isinstance(r, (list, tuple)) else r)(G(v))
  return ca.horzcat(*[(f(v0 + h*e) - f(v0 - h*e))/(2*h) for e in ca.horzsplit(ca.DM.eye(v0.numel()))])

print("Will test these integrators:")
for cl, t, options in integrators:
  print(cl, " : ", t)

class Integrationtests(casadiTestCase):

  @slow()
  @memory_heavy()
  def test_full(self):
    num = self.num
    tc = ca.DM(n.linspace(0,num['tend'],100))

    t=ca.SX.sym("t")
    q=ca.SX.sym("q")
    p=ca.SX.sym("p")

    dae = {'t':t, 'x':q, 'p':p, 'ode':q/p*t**2}
    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    opts["fsens_err_con"] = True
    #opts["verbose"] = True
    integrator = ca.integrator("integrator", "cvodes", dae, 0, 2.3, opts)
    with self.assertInException("No stats available"):
        integrator.stats()
    tf = 2.3

    solution = ca.Function("solution", {'x0':q, 'p':p, 'xf':q*ca.exp(tf**3/(3*p))},
                        ca.integrator_in(), ca.integrator_out())
    f_in = {}
    f_in["x0"]=0.3
    f_in["p"]=0.7

    self.checkfunction(integrator,solution,inputs=f_in,digits=6)

  def test_tools_trivial(self):
    num = self.num

    x = ca.SX.sym("x")
    p = ca.SX.sym("p",0)
    tf = ca.SX.sym("tf")
    f=ca.Function("f", [x,p], [x])

    for integrator in [
         ca.simpleRK(f),
         ca.simpleIRK(f),
         ca.simpleIntegrator(f)
       ]:

      solution = ca.Function("solution", [x,p,tf], [x*ca.exp(tf)], ["x0","p","h"], ["xf"])

      f_in = {"x0":1,"h":1}
      self.checkfunction(integrator,solution,inputs=f_in,digits=3)

  @slow()
  @memory_heavy()
  def test_tools(self):
    num = self.num

    t=ca.SX.sym("t")
    q=ca.SX.sym("q")
    p=ca.SX.sym("p")
    t0=ca.SX.sym("t0")
    q0=ca.SX.sym("q0")
    f=ca.Function("f", [ca.vertcat(*(q,t)),p],[ca.vertcat(*(q/p*t**2,1))])
    for integrator in [
            ca.simpleRK(f,500),
            ca.simpleIRK(f,50),
            ca.simpleIntegrator(f)
            ]:

      solution = ca.Function('solver', {"x0": ca.vertcat(*(q0,t0)),"p":p,"h":t, "xf":ca.vertcat(*[q0*ca.exp(((t0+t)**3-t0**3)/(3*p)),t0+t])}, integrator.name_in(), integrator.name_out())

      f_in = {}
      f_in["x0"]=ca.DM([0.3,0])
      f_in["p"]=0.7
      f_in["h"]=1

      self.checkfunction(integrator,solution,inputs=f_in,digits=3)

  @memory_heavy()
  def test_jac(self):
    self.message("Test exact jacobian #536")
    # This test is not automized, but works by inspection only.
    # To activate, recompile after ucnommenting the printout lines in cvodes.c, near "Used for validating casadi#536"
    #return
    ca.DM.set_precision(18)

    tstart = ca.SX.sym("tstart")
    tend = ca.SX.sym("tend")

    integrators = [
              ("idas",["dae","ode"],{"abstol": 1e-9,"reltol":1e-9,"fsens_err_con": True,"calc_ic":True,"calc_icB":True}),
              ("cvodes",["ode"],{"abstol": 1e-5,"reltol":1e-5,"fsens_err_con": False,"quad_err_con": False})
              ]

    def variations(p_features, din, dout, *args):
      if "ode" in p_features:
        p_features_ = copy.copy(p_features)
        p_features_[p_features.index("ode")] = "dae"
        din_ = copy.copy(din)
        dout_ = copy.copy(dout)
        z = ca.SX.sym("x", din_["x"].shape)
        din_["z"] = z
        dout_["ode"] = z
        dout_["alg"] = ( dout["ode"] - z) * (-0.8)

        yield (p_features, din, dout) + tuple(args)
        yield (p_features_, din_, dout_) + tuple(args)
      else:
        yield (p_features, din, dout) + tuple(args)

    def checks():
      Ns = 1

      x  = ca.SX.sym("x")
      t = ca.SX.sym("t")

      ti = (0,0.9995)
      pointA = {'x0': 1}

      si = {'x0':x}

      sol = {'xf': x*ca.exp(tend)}

      yield (["ode"],{'x':x,'t':t},{'ode':x},si,sol,pointA,ti)


    refXF = refRXF = None

    for tt in checks():
      for p_features, din, dout, solutionin, solution, point, (tstart_, tend_) in variations(*tt):
        for Integrator, features, options in integrators:
          self.message(Integrator)
          x = ca.SX.sym("x")
          dummyIntegrator = ca.integrator("dummyIntegrator", Integrator, {'x':x, 'ode':x})
          if p_features[0] in features:

            dae = din.copy()
            dae.update(dout)

            for k in list(solution.keys()):
              solution[k] = ca.substitute(solution[k],ca.vertcat(*[tstart,tend]),ca.vertcat(*[tstart_,tend_]))

            sol = solutionin.copy()
            sol.update(solution)
            fs = ca.Function("fs", sol, ca.integrator_in(), ca.integrator_out())

            def itoptions():
                yield {"newton_scheme": "direct"}
                yield {"newton_scheme": "gmres"}
                yield {"newton_scheme": "bcgstab"}
                # yield {"newton_scheme": "tfqmr"} # Bug in Sundials?

            for f_options in itoptions():
                message = "f_options: %s" % str(f_options)
                print(message)
                opts = {}
                for op in (options, f_options):
                  for (k,v) in list(op.items()):
                    opts[k] = v
                integrator = ca.integrator("integrator", Integrator, dae, tstart_, tend_, opts)
                integrator_in = {}
                for k,v in list(point.items()):
                  if not fs.sparsity_in(k).is_empty():
                    integrator_in[k]=v

                integrator_out = integrator(**integrator_in)
                fs_out = fs(**integrator_in)
                print("res=",integrator_out["xf"]-fs_out["xf"], fs_out["xf"])
                stats = integrator.stats()

                print(stats)

                self.assertTrue(stats["nstepsB"]<2500)


                self.assertTrue(stats["nsteps"]<1500)
                self.assertTrue(stats["nstepsB"]<2500)
                self.assertTrue(stats["nlinsetups"]<100)

                self.assertTrue(stats["nlinsetupsB"]<250)

    ca.DM.set_precision(6)

  @memory_heavy()
  def test_lsolvers(self):
    self.message("Test different linear solvers")

    tstart = ca.SX.sym("tstart")
    tend = ca.SX.sym("tend")

    integrators = [
              ("idas",["dae","ode"],{"abstol": 1e-9,"reltol":1e-9,"fsens_err_con": True,"calc_ic":True,"calc_icB":True}),
              ("cvodes",["ode"],{"abstol": 1e-15,"reltol":1e-15,"fsens_err_con": True,"quad_err_con": False})
              ]

    def checks():
      t=ca.SX.sym("t")
      x=ca.SX.sym("x")
      p=ca.SX.sym("p")
      dp=ca.SX.sym("dp")

      z=ca.SX.sym("z")
      solutionin = {'x0':x, 'p': p}
      pointA = {'x0':7.1,'p': 2}
      ti = (0.2,2.3)
      yield (["dae"],{'x': x, 'z': z},{'alg': x-z, 'ode': z},solutionin,{'xf':x*ca.exp(tend-tstart)},pointA,ti)
      if not(args.run_slow): return
      yield (["dae"],{'x': x, 'z': z},{'alg': x-z, 'ode': z},solutionin,{'xf':x*ca.exp(tend-tstart)},pointA,ti)
      yield (["ode"],{'x': x},{'ode': x},solutionin,{'xf':x*ca.exp(tend-tstart)},pointA,ti)

      A=array([1,0.1])
      p0 = 1.13

      q=ca.SX.sym("y",2,1)
      y0=q[0]
      yc0=dy0=q[1]
      p=ca.SX.sym("p",1,1)

      s1=(2*y0-ca.log(yc0**2/p+1))/2-ca.log(ca.cos(ca.arctan(yc0/ca.sqrt(p))+ca.sqrt(p)*(tend-tstart)))
      s2=ca.sqrt(p)*ca.tan(ca.arctan(yc0/ca.sqrt(p))+ca.sqrt(p)*(tend-tstart))
      yield (["ode"],{'x':q,'p':p},{'ode': ca.vertcat(*[q[1],p[0]+q[1]**2 ])},{'x0':q, 'p': p} ,{'xf': ca.vertcat(*[s1,s2])},{'x0': A, 'p': p0},(0,0.4) )

    for p_features, din, dout, solutionin, solution, point, (tstart_, tend_) in checks():

      for Integrator, features, options in integrators:
        self.message(Integrator)
        x = ca.SX.sym("x")
        dummyIntegrator = ca.integrator("dummyIntegrator", Integrator, {'x':x, 'ode':x})
        if p_features[0] in features:

          dae = din.copy()
          dae.update(dout)

          for k in list(solution.keys()):
            solution[k] = ca.substitute(solution[k],ca.vertcat(*[tstart,tend]),ca.vertcat(*[tstart_,tend_]))

          sol = solutionin.copy()
          sol.update(solution)
          fs = ca.Function("fs", sol, ca.integrator_in(), ca.integrator_out())

          def itoptions():
            yield {"newton_scheme": "direct"}
            yield {"newton_scheme": "gmres"}
            yield {"newton_scheme": "bcgstab"}
            #yield {"newton_scheme": "tfqmr"} # Bug in Sundials?

          for f_options in itoptions():
              message = "f_options: %s" % str(f_options)
              print(message)

              opts = {}
              for op in (options, f_options):
                 for (k,v) in list(op.items()):
                    opts[k] = v
              integrator = ca.integrator("integrator", Integrator, dae, tstart_, tend_, opts)
              integrator_in = {}

              for k,v in list(point.items()):
                if not integrator.sparsity_in(k).is_empty():
                  integrator_in[k]=v

              self.checkfunction(integrator,fs,inputs=integrator_in,gradient=False,hessian=False,sens_der=False,evals=False,digits=4,digits_sens=4,failmessage=message,verbose=False)
              self.check_pure(integrator,inputs=integrator_in)
              self.check_serialize(integrator,inputs=integrator_in)

  @memory_heavy()
  def test_X(self):
    self.message("Extensive integrator tests")

    num=self.num
    tstart = ca.SX.sym("tstart")
    tend = ca.SX.sym("tend")


    for Integrator, features, options in integrators:
      self.message(Integrator)


      def variations(p_features, din, dout, *args):
        if "ode" in p_features:
          p_features_ = copy.copy(p_features)
          p_features_[p_features.index("ode")] = "dae"
          din_ = copy.copy(din)
          dout_ = copy.copy(dout)
          z = ca.SX.sym("z", din_["x"].shape)
          din_["z"] = z
          dout_["ode"] = z
          dout_["alg"] = ( dout["ode"] - z) * (-0.8)

          yield (p_features, din, dout) + tuple(args)
          yield (p_features_, din_, dout_) + tuple(args)
        else:
          yield (p_features, din, dout) + tuple(args)

      def checks():
        x0=num['q0']
        p_=num['p']
        t=ca.SX.sym("t")
        x=ca.SX.sym("x")
        p=ca.SX.sym("p")
        dp=ca.SX.sym("dp")

        z=ca.SX.sym("z")

        si = {'x0':x, 'p': p}
        pointA = {'x0':x0,'p': p_}

        ti = (0.2,num['tend'])
        yield (["ode"],{'x':x},{'ode': 0},si,{'xf':x},pointA,ti)
        yield (["ode"],{'x':x},{'ode': 1},si,{'xf':x+(tend-tstart)},pointA,ti)
        yield (["ode"],{'x':x},{'ode': x},si,{'xf':x*ca.exp(tend-tstart)},pointA,ti)
        yield (["ode"],{'x':x,'t':t},{'ode': t},si,{'xf':x+(tend**2/2-tstart**2/2)},pointA,ti)
        yield (["ode"],{'x':x,'t':t},{'ode': x*t},si,{'xf':x*ca.exp(tend**2/2-tstart**2/2)},pointA,ti)
        yield (["ode"],{'x':x,'p':p},{'ode': x/p},si,{'xf':x*ca.exp((tend-tstart)/p)},pointA,ti)
        if not(args.run_slow): return
        yield (["ode"],{'x':x},{'ode': x,'quad':0},si,{'qf':0},pointA,ti)
        yield (["ode"],{'x':x},{'ode': x,'quad':1},si,{'qf':(tend-tstart)},pointA,ti)
        yield (["ode"],{'x':x},{'ode': 0,'quad':x},si,{'qf':x*(tend-tstart)},pointA,ti)
        yield (["ode"],{'x':x},{'ode': x,'quad':x},si,{'qf':x*(ca.exp(tend-tstart)-1)},pointA,ti)
        yield (["ode"],{'x':x,'t':t},{'ode': x,'quad':t},si,{'qf':(tend**2/2-tstart**2/2)},pointA,ti)
        yield (["ode"],{'x':x,'t':t},{'ode': x,'quad':x*t},si,{'qf':x*(ca.exp(tend-tstart)*(tend-1)-(tstart-1))},pointA,ti)
        yield (["ode"],{'x':x,'p':p},{'ode': x,'quad':x/p},si,{'qf':x*(ca.exp((tend-tstart))-1)/p},pointA,ti)
        yield (["dae"],{'x':x,'z':z},{'ode':z,'alg': -0.8*(z-x),'quad': z},si,{'qf':x*(ca.exp(tend-tstart)-1)},pointA,ti)


        A=array([1,0.1])
        p0 = 1.13

        q=ca.SX.sym("y",2,1)
        y0=q[0]
        yc0=dy0=q[1]
        p=ca.SX.sym("p",1,1)

        s1=(2*y0-ca.log(yc0**2/p+1))/2-ca.log(ca.cos(ca.arctan(yc0/ca.sqrt(p))+ca.sqrt(p)*(tend-tstart)))
        s2=ca.sqrt(p)*ca.tan(ca.arctan(yc0/ca.sqrt(p))+ca.sqrt(p)*(tend-tstart))
        yield (["ode"],{'x':q,'p':p},{'ode': ca.vertcat(*[q[1],p[0]+q[1]**2 ])},{'x0':q, 'p': p} ,{'xf': ca.vertcat(*[s1,s2])},{'x0': A, 'p': p0},(0,0.4) )

      for tt in checks():
        print(tt)
        for p_features, din, dout, solutionin, solution, point, (tstart_, tend_) in variations(*tt):
          if p_features[0] in features:
            message = "%s: %s => %s, explicit (%s) tstart = %f" % (Integrator,str(din),str(dout),str(solution),tstart_)
            print(message)

            dae = din.copy()
            dae.update(dout)

            for k in list(solution.keys()):
              solution[k] = ca.substitute(solution[k],ca.vertcat(*[tstart,tend]),ca.vertcat(*[tstart_,tend_]))

            sol = solutionin.copy()
            sol.update(solution)
            fs = ca.Function("fs", sol, ca.integrator_in(), ca.integrator_out())

            opts = dict(options)
            if Integrator in ('cvodes', 'idas'):
              opts["abstol"] = 1e-9
              opts["reltol"] = 1e-9
            if Integrator=='idas':
              opts["init_xdot"] = ca.DM(point["x0"]).nonzeros()
              opts["calc_icB"] = True
              opts["augmented_options"] = {"init_xdot":None, "abstol":1e-9,"reltol":1e-9}
            integrator = ca.integrator("integrator", Integrator, dae, tstart_, tend_, opts)
            integrator_in = {}

            for k,v in list(point.items()):
              if not integrator.sparsity_in(k).is_empty():
                integrator_in[k]=v

            self.checkfunction(integrator,fs,inputs=integrator_in,gradient=False,hessian=False,sens_der=False,evals=False,digits=4,digits_sens=4,failmessage=message,verbose=False)
            

            if "kinsol" not in str(opts):
              print(integrator)
              self.check_serialize(integrator,inputs=integrator_in)
              
  @memory_heavy()
  def test_multipoint(self):
  
    for Integrator, features, options in integrators:
      self.message(Integrator)
      
      if str(Integrator) in ["cvodes","idas"]: continue
      
      x = ca.MX.sym("x")
      p = ca.MX.sym("p")
      z = ca.MX.sym("z")
      u = ca.MX.sym("u")
      
      if "dae" in features:
          dae = { 'x':x, 'z':z, 'p':p, 'u': u, 'ode':z+p*u, 'alg':z*ca.cos(z)-x, 'quad': x**2}
      else:
          dae = { 'x':x, 'p':p, 'u': u, 'ode':ca.cos(x+p*u), 'quad': x**2}
          
      integratorA = ca.integrator("integrator", Integrator, dae, 0, 0.1, options)
      integratorB = ca.integrator("integrator", Integrator, dae, 0.1, 0.2, options)
   
      integratorI = ca.integrator("integrator", Integrator, dae, 0, [0.1,0.2], options)
      
      arg_test = {"x0":0,"p":0.15,"u":ca.DM([1,1.1]).T}
    
      args = integratorI.convert_in(integratorI.mx_in())
      u_split = ca.horzsplit_n(args["u"],2)
      
      args1 = dict(args)
      args2 = dict(args)
      args1["u"] = u_split[0]
      args2["u"] = u_split[1]
      res1 = integratorA(**args1)
      args2["x0"] = res1["xf"]
      res2 = integratorB(**args2)
      res2["qf"] = ca.horzcat(res1["qf"],res1["qf"]+res2["qf"])
      res2["xf"] = ca.horzcat(res1["xf"],res2["xf"])
      res2["zf"] = ca.horzcat(res1["zf"],res2["zf"])
      
      res2.update(args)
      I = ca.Function("I",res2,integratorI.name_in(),integratorI.name_out())

      I(**arg_test)      
      integratorI(**arg_test)
      
      self.checkfunction(integratorI,I,inputs={"x0":0,"p":0.15,"u":ca.DM([1,1.1]).T})
            

  def setUp(self):
    # Reference solution is x0 e^((t^3-t0^3)/(3 p))
    t=ca.SX.sym("t")
    x=ca.SX.sym("x")
    p=ca.SX.sym("p")
    dae={'t':t, 'x':x, 'p':p, 'ode':x/p*t**2}
    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    opts["fsens_err_con"] = True
    #opts["verbose"] = True
    integrator = ca.integrator("integrator", "cvodes", dae, 0, 2.3, opts)
    q0   = ca.MX.sym("q0")
    par  = ca.MX.sym("p")

    qend = integrator(x0=q0, p=par)["xf"]

    qe=ca.Function("qe", [q0,par],[qend])
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
    qe=ca.Function("qe", [q0,par],[qend[0]])

    f = ca.Function("f", [q0], [qe(q0,ca.MX(num['p']))])
    f_out = f([num['q0']])
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(f_out[0][0],q0*ca.exp(tend**3/(3*p)),9,"Evaluation output mismatch")

  def test_issue92b(self):
    self.message("regression check for issue 92")
    t=ca.SX.sym("t")
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")
    dae = {'t':t, 'x':ca.vertcat(*[x,y]), 'ode':ca.vertcat(*[x,(1+1e-9)*x])}
    opts = {}
    opts["fsens_err_con"] = True
    integrator = ca.integrator("integrator", "cvodes", dae, opts)
    integrator_in = {}

    # Pass inputs
    integrator_in["x0"]=ca.DM([1,0])
    ## Integrate
    integrator_out = integrator(**integrator_in)
    # print result
    print(integrator_out["xf"])

  def test_issue92(self):
    self.message("regression check for issue 92")
    t=ca.SX.sym("t")
    x=ca.SX.sym("x")
    var = ca.MX.sym("var",2,1)

    q = ca.vertcat(*[x,ca.SX.sym("problem")])

    dq=ca.vertcat(*[x,x])
    dae = {'t':t, 'x':q, 'ode':dq}
    opts = {}
    opts["fsens_err_con"] = True
    opts["reltol"] = 1e-12
    integrator = ca.integrator("integrator", "cvodes", dae, opts)

    qend = integrator(x0=var)["xf"]

    f = ca.Function("f", [var],[qend[0]])

    J = jacobian_old(f, 0, 0)
    J_out = J([1,0])
    print("jac=",J_out[0].nz[0]-ca.exp(1))
    self.assertAlmostEqual(J_out[0][0,0],ca.exp(1),5,"Evaluation output mismatch")

  def test_eval(self):
    self.message('CVodes integration: evaluation')
    num=self.num
    qe=self.qe
    qe_out = qe([num['q0']], [num['p']])
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(qe_out[0][0],q0*ca.exp(tend**3/(3*p)),9,"Evaluation output mismatch")

  def test_jac1(self):
    self.message('CVodes integration: jacobian to q0')
    num=self.num
    J = jacobian_old(self.qe, 0, 0)
    J_out = J([num['q0']], [num['p']])
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(J_out[0][0],ca.exp(tend**3/(3*p)),9,"Evaluation output mismatch")

  def test_jac2(self):
    self.message('CVodes integration: jacobian to p')
    num=self.num
    J = jacobian_old(self.qe, 1, 0)
    J_out = J([num['q0']], [num['p']])
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(J_out[0][0],-(q0*tend**3*ca.exp(tend**3/(3*p)))/(3*p**2),9,"Evaluation output mismatch")

  def test_bug_repeat(self):
    num={'tend':2.3,'q0':[0,7.1,7.1],'p':2}
    self.message("Bug that appears when rhs contains repeats")
    A=array([1,0.1,1])
    p0 = 1.13
    y0=A[0]
    yc0=dy0=A[1]
    te=0.4

    t=ca.SX.sym("t")
    q=ca.SX.sym("y",3,1)
    p=ca.SX.sym("p")

    dh = p+q[0]**2
    dae = {'x':q, 'p':p, 't':t, 'ode':ca.vertcat(*[dh ,q[0],dh])}

    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    #opts["verbose"] = True
    opts["fsens_err_con"] = True
    opts["steps_per_checkpoint"] = 10000
    integrator = ca.integrator("integrator", "cvodes", dae, 0, te, opts)

    q0   = ca.MX.sym("q0",3,1)
    par  = ca.MX.sym("p",1,1)
    qend = integrator(x0=q0, p=par)["xf"]
    qe=ca.Function("qe", [q0,par],[qend])

    #J=self.qe.jacobian_old(2, 0)
    J = jacobian_old(qe, 0, 0)
    J_out = J(A, p0)
    outA=J_out[0].full()
    dae={'x':q, 'p':p, 't':t, 'ode':ca.vertcat(*[dh ,q[0],(1+1e-9)*dh])}

    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    #opts["verbose"] = True
    opts["fsens_err_con"] = True
    opts["steps_per_checkpoint"] = 10000
    integrator = ca.integrator("integrator", "cvodes", dae, 0, te, opts)

    q0   = ca.MX.sym("q0",3,1)
    par  = ca.MX.sym("p",1,1)
    qend = integrator(x0=q0, p=par)["xf"]
    qe=ca.Function("qe", [q0,par],[qend])

    #J=self.qe.jacobian_old(2)
    J = jacobian_old(qe, 0, 0)
    J_out = J(A, p0)
    outB=J_out[0].full()
    print(outA-outB)

  def test_hess3(self):
    self.message('CVodes integration: hessian to p: Jacobian of integrator.jacobian')
    num=self.num
    J = jacobian_old(self.integrator, self.integrator.index_in("p"),self.integrator.index_out("xf"))
    H = jacobian_old(J, J.index_in("p"), 0)
    H_in = {}
    H_in["x0"]=num['q0']
    H_in["p"]=num['p']
    H_out = H(**H_in)
    num=self.num
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(H_out["jac_jac_xf_p_p"][0],(q0*tend**6*ca.exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*ca.exp(tend**3/(3*p)))/(3*p**3),9,"Evaluation output mismatch")

  def test_hess4(self):
    self.message('CVodes integration: hessian to p: Jacobian of integrator.jacobian indirect')
    num=self.num
    J = jacobian_old(self.integrator, self.integrator.index_in("p"),self.integrator.index_out("xf"))

    q0=ca.MX.sym("q0")
    p=ca.MX.sym("p")
    Ji = J(x0=q0,p=p)
    Ji["q0"] = q0
    Ji["p"] = p
    Ji = ca.Function("Ji", Ji, ["q0", "p"], J.name_out())
    H = jacobian_old(Ji, 1, 0)
    H_out = H([num['q0']], [num['p']])
    num=self.num
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(H_out[0][0],(q0*tend**6*ca.exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*ca.exp(tend**3/(3*p)))/(3*p**3),9,"Evaluation output mismatch")

  def test_hess5(self):
    self.message('CVodes integration: hessian to p in an MX tree')
    num=self.num
    q0=ca.MX.sym("q0")
    p=ca.MX.sym("p")

    sol = self.integrator(x0=q0,p=p)
    sol["q0"] = q0
    sol["p"] = p
    JT = ca.Function("JT", [q0,p],[ca.jacobian(sol['xf'], sol['p']).T])
    JT_out = JT([num['q0']], [num['p']])
    print(JT_out)

    H = jacobian_old(JT, 1, 0)
    H_out = H([num['q0']], [num['p']])
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(H_out[0],(q0*tend**6*ca.exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*ca.exp(tend**3/(3*p)))/(3*p**3),9,"Evaluation output mismatch")

  def test_hess6(self):
    self.message('CVodes integration: hessian to p in an MX tree')
    num=self.num
    q0=ca.MX.sym("q0")
    p=ca.MX.sym("p")

    sol = self.integrator(x0=q0,p=p)
    sol["q0"] = q0
    sol["p"] = p
    qe = ca.Function("qe", sol, ["q0", "p"], ca.integrator_out())

    H = hessian_old(qe, 1, 0)
    H_out = H([num['q0']], [num['p']])
    num=self.num
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(H_out[0][0],(q0*tend**6*ca.exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*ca.exp(tend**3/(3*p)))/(3*p**3),7,"Evaluation output mismatch")

  def test_glibcbug(self):
    self.message("former glibc error")
    A=array([2.3,4.3,7.6])
    B=array([[1,2.3,4],[-2,1.3,4.7],[-2,6,9]])

    te=0.7
    t=ca.SX.sym("t")
    q=ca.SX.sym("q",3,1)
    p=ca.SX.sym("p",9,1)
    dae = {'t':t, 'x':q, 'p':p, 'ode':c.reshape(p,3,3) @ q}
    opts = {}
    opts["fsens_err_con"] = True
    opts["steps_per_checkpoint"] = 1000
    integrator = ca.integrator("integrator", "cvodes", dae, 0, te, opts)
    q0   = ca.MX.sym("q0",3,1)
    par  = ca.MX.sym("p",9,1)
    qend = integrator(x0=q0, p=par)["xf"]
    qe = jacobian_old(integrator, integrator.index_in("p"),integrator.index_out("xf"))
    qe = qe(x0=q0,p=par)['jac_xf_p']
    qef=ca.Function("qef", [q0,par],[qe])
    qef_out = qef(A, B.ravel())

  def test_linear_system(self):
    self.message("Linear ODE")
    if not(scipy_available):
        return
    A=array([2.3,4.3,7.6])
    B=array([[1,2.3,4],[-2,1.3,4.7],[-2,6,9]])
    te=0.7
    Be=expm(B*te)
    t=ca.SX.sym("t")
    q=ca.SX.sym("q",3,1)
    p=ca.SX.sym("p",9,1)

    dae = {'t':t, 'x':q, 'p':p, 'ode':c.reshape(p,3,3) @ q}
    opts = {}
    opts["fsens_err_con"] = True
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    #opts["verbose"] = True
    opts["steps_per_checkpoint"] = 10000
    integrator = ca.integrator("integrator", "cvodes", dae, 0, te, opts)

    q0   = ca.MX.sym("q0",3,1)
    par  = ca.MX.sym("p",9,1)
    qend = integrator(x0=q0, p=par)["xf"]
    qe=ca.Function("qe", [q0,par],[qend])
    qendJ = jacobian_old(integrator, integrator.index_in("x0"),integrator.index_out("xf"))
    qendJ = qendJ(x0=q0,p=par)['jac_xf_x0']

    qeJ=ca.Function("qeJ", [q0,par],[qendJ])

    qendJ2 = jacobian_old(integrator, integrator.index_in("x0"),integrator.index_out("xf"))
    qendJ2 = qendJ2(x0=q0,p=par)['jac_xf_x0']

    qeJ2=ca.Function("qeJ2", [q0,par],[qendJ2])
    qe_out = qe(A, ca.vec(B))
    self.checkarray(numpy.dot(Be,A)/1e3,qe_out/1e3,"jacobian('x0','xf')")
    qeJ_out = qeJ(A, ca.vec(B))
    self.checkarray(qeJ_out/1e3,Be/1e3,"jacobian('x0','xf')")

    qeJ2_out = qeJ2(A, ca.vec(B))

    return # this should return identical zero
    H = jacobian_old(qeJ, 0, 0)
    H_out = H(A, ca.vec(B))
    print(array(H_out[0]))

  def test_mathieu_system(self):
    self.message("Mathieu ODE")
    A=array([0.3,1.2])
    B=array([1.3,4.3,2.7])
    te=0.7

    t=ca.SX.sym("t")
    q=ca.SX.sym("q",2,1)
    p=ca.SX.sym("p",3,1)

    dae = {'x':q, 'p':p, 't':t, 'ode':ca.vertcat(*[q[1],(p[0]-2*p[1]*ca.cos(2*p[2]))*q[0]])}
    opts = {}
    opts["fsens_err_con"] = True
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    #opts["verbose"] = True
    opts["steps_per_checkpoint"] = 10000
    integrator = ca.integrator("integrator", "cvodes", dae, 0, te, opts)

    q0   = ca.MX.sym("q0",2,1)
    par  = ca.MX.sym("p",3,1)
    qend = integrator(x0=q0, p=par)["xf"]
    qe=ca.Function("qe", [q0,par],[qend])
    qendJ = jacobian_old(integrator, integrator.index_in("x0"), integrator.index_out("xf"))
    qendJ =qendJ(x0=q0,p=par)['jac_xf_x0']
    qeJ=ca.Function("qeJ", [q0,par],[qendJ])
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

    t=ca.SX.sym("t")
    q=ca.SX.sym("y",2,1)
    p=ca.SX.sym("p",1,1)
    # y
    # y'
    dae = {'x':q, 'p':p, 't':t, 'ode':ca.vertcat(*[q[1],p[0]+q[1]**2])}
    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    #opts["verbose"] = True
    opts["steps_per_checkpoint"] = 10000
    opts["fsens_err_con"] = True
    integrator = ca.integrator("integrator", "cvodes", dae, 0, te, opts)

    t0   = ca.MX(0)
    tend = ca.MX(te)
    q0   = ca.MX.sym("q0",2,1)
    par  = ca.MX.sym("p",1,1)
    qend = integrator(x0=q0, p=par)["xf"]
    qe=ca.Function("qe", [q0,par],[qend])
    qendJ = jacobian_old(integrator, integrator.index_in("x0"), integrator.index_out("xf"))
    qendJ = qendJ(x0=q0, p=par)['jac_xf_x0']
    qeJ=ca.Function("qeJ", [q0,par],[qendJ])
    qe_out = qe(A, p0)

    print(qe_out[0])
    print(qe_out[1])

    self.assertAlmostEqual(qe_out[0],(2*y0-ca.log(yc0**2/p0+1))/2-ca.log(ca.cos(ca.arctan(yc0/ca.sqrt(p0))+ca.sqrt(p0)*te)),11,"Nonlin ODE")
    self.assertAlmostEqual(qe_out[1],ca.sqrt(p0)*ca.tan(ca.arctan(yc0/ca.sqrt(p0))+ca.sqrt(p0)*te),11,"Nonlin ODE")

    qeJ_out = qeJ(A, p0)

    Jr = array([[1,(ca.sqrt(p0)*ca.tan(ca.sqrt(p0)*te+ca.arctan(dy0/ca.sqrt(p0)))-dy0)/(dy0**2+p0)],[0,(p0*ca.tan(ca.sqrt(p0)*te+ca.arctan(dy0/ca.sqrt(p0)))**2+p0)/(dy0**2+p0)]])
    self.checkarray(qeJ_out,Jr,"jacobian of Nonlin ODE")

    Jf = jacobian_old(qe, 0,0)
    Jf_out = Jf(A, p0)
    self.checkarray(Jf_out[0],Jr,"Jacobian of Nonlin ODE")

    Jf = jacobian_old(qe, 0, 0)
    Jf_out = Jf(A, p0)
    self.checkarray(Jf_out[0],Jr,"Jacobian of Nonlin ODE")

    Jr = numpy.array([[(ca.sqrt(p0)*(te*yc0**2-yc0+p0*te)*ca.tan(ca.arctan(yc0/ca.sqrt(p0))+ca.sqrt(p0)*te)+yc0**2)/(2*p0*yc0**2+2*p0**2)],[(ca.sqrt(p0)*((te*yc0**2-yc0+p0*te)*ca.tan(ca.arctan(yc0/ca.sqrt(p0))+ca.sqrt(p0)*te)**2+te*yc0**2-yc0+p0*te)+(yc0**2+p0)*ca.tan(ca.arctan(yc0/ca.sqrt(p0))+ca.sqrt(p0)*te))/(ca.sqrt(p0)*(2*yc0**2+2*p0))]])

    Jf = jacobian_old(qe, 1, 0)
    Jf_out = Jf(A, p0)
    self.checkarray(Jf_out[0],Jr,"Jacobian of Nonlin ODE")
    Jf = jacobian_old(qe, 1, 0)
    Jf_out = Jf(A, p0)
    self.checkarray(Jf_out[0],Jr,"Jacobian of Nonlin ODE")

    qendJ = jacobian_old(integrator, integrator.index_in("p"),integrator.index_out("xf"))
    qendJ = qendJ(x0=q0,p=par)['jac_xf_p']
    qeJ=ca.Function("qeJ", [q0,par],[qendJ])

    qeJ_out = qeJ(A, p0)

    self.checkarray(qeJ_out,Jr,"jacobian of Nonlin ODE")

    self.check_thread_safety(qeJ,inputs=[A, p0])

    qeJf=ca.Function("qeJf", [q0,par],[ca.vec(qeJ(q0,par))])

    H = jacobian_old(qeJf, 0, 0)
    H_out = H(A, p0)
    def sec(x):
      return 1.0/ca.cos(x)
    Hr = array([[0,0],[0,-(2*yc0*ca.tan(ca.arctan(yc0)+te))/(yc0**4+2*yc0**2+1)+sec(ca.arctan(yc0)+te)**2/(yc0**4+2*yc0**2+1)+(2*yc0**2)/(yc0**4+2*yc0**2+1)-1/(yc0**2+1)],[0,0],[0,-(2*yc0*ca.tan(ca.arctan(yc0)+te)**2)/(yc0**4+2*yc0**2+1)+(2*sec(ca.arctan(yc0)+te)**2*ca.tan(ca.arctan(yc0)+te))/(yc0**4+2*yc0**2+1)-(2*yc0)/(yc0**4+2*yc0**2+1)]])
    print(array(H_out[0]))
    print(Hr)

  def test_missing_symbols(self):
    x = ca.SX.sym('x'); z = ca.SX.sym('z'); p = ca.SX.sym('p')
    dae = {'x':x, 'z':z, 'ode':z+p, 'alg':z*ca.cos(z)-x} # p forgotten here
    with self.assertInException("[p] are free"):
      ca.integrator('F', 'idas', dae)

  def test_hessian2D(self):
    self.message("hessian")
    N=2

    x0_ = ca.DM([1,0.1])
    A_  = ca.DM([[3,1],[0.74,4]])

    A = ca.SX.sym("A",N,N)
    x = ca.SX.sym("x",N)

    dae = {'x':x, 'p':ca.vec(A), 'ode':A @ x}
    I = ca.integrator("I", "cvodes", dae, {"fsens_err_con": True, 'reltol' : 1e-12})
    I_in =  {}
    I_in["x0"]=x0_
    I_in["p"]=ca.vec(A_)
    I_out = I(**I_in)

    q0=ca.MX.sym("q0",N)
    p=ca.MX.sym("p",N*N)

    sol = I(x0=q0,p=p)
    sol["q0"] = q0
    sol["p"] = p
    qe = ca.Function("qe", sol, ["q0", "p"], ca.integrator_out())

    JT = ca.Function("JT", [q0,p],[ca.jacobian(sol['xf'], sol['p']).T])

    H = jacobian_old(JT, 1, 0)
    H_out = H(x0_, ca.vec(A_))


    H1 = H_out[0]

  def test_collocationPoints(self):
    self.message("collocation points")
    with self.assertRaises(Exception):
      ca.collocation_points(0,"radau")
    with self.assertRaises(Exception):
      ca.collocation_points(10,"radau")
    with self.assertRaises(Exception):
      ca.collocation_points(0,"legendre")
    with self.assertRaises(Exception):
      ca.collocation_points(10,"legendre")
    with self.assertRaises(Exception):
      ca.collocation_points(1,"foo")

    for k in range(1,10):
      r = [0] + ca.collocation_points(k,"radau")
      self.assertEqual(len(r),k+1)
      self.checkarray(ca.DM(r[-1]),ca.DM([1]))
    for k in range(1,10):
      r = [0] + ca.collocation_points(k,"legendre")
      self.assertEqual(len(r),k+1)

  @memory_heavy()
  def test_thread_safety(self):
    x = ca.MX.sym('x')
    t = ca.MX.sym('t')
    p = ca.MX.sym('p')
    for Integrator, features, options in integrators:
      intg = ca.integrator("Integrator",Integrator,{"x":x,"ode":-x}, options)

      intg_par = intg.map(40, 'thread',4)
      res = intg_par(x0=numpy.linspace(0, 10, 40))
      self.checkarray(ca.norm_inf(res["xf"].T-ca.exp(-1)*numpy.linspace(0, 10, 40)),0, digits=5)


      intg = ca.integrator("Integrator",Integrator,{"x":x,"ode":-x}, options)

      intg_par = intg.map(40, 'thread',2)
      res = intg_par(x0=numpy.linspace(0, 10, 40))
      self.checkarray(ca.norm_inf(res["xf"].T-ca.exp(-1)*numpy.linspace(0, 10, 40)),0, digits=5)

  def test_simplify_zdim(self):
    x = ca.MX.sym("x")
    intg = ca.integrator("intg","rk",{"x":x,"ode":x**2},{"simplify":True})

    self.assertTrue(intg.nnz_out("zf")==0)

  @requires_integrator('cvodes')
  def test_step_options_cvodes(self):
    x = ca.SX.sym("x")
    opts = {
      "step0":    1e-4,
      "min_step_size": 1e-4,
      "max_step_size": 1.1e-4
    }

    I = ca.integrator("I","cvodes",{"x":x,"ode":ca.sin((10*x)**2)}, 0.0, 0.5, opts)
    I(x0=0)
    stats = I.stats()

    self.assertTrue(int(0.5/1e-4)>=stats["nsteps"]>=int(0.5/1.1e-4))

  @requires_integrator('idas')
  def test_step_options_idas(self):
    x = ca.SX.sym("x")
    opts = {
      "step0":    1e-4,
      "max_step_size": 1.1e-4
    }

    I = ca.integrator("I","cvodes",{"x":x,"ode":ca.sin((10*x)**2)}, 0.0, 0.5, opts)
    I(x0=0)
    stats = I.stats()

    print(stats["nsteps"])
    self.assertTrue(stats["nsteps"]>=int(0.5/1.1e-4))

  def test_adjoint_without_parameters(self):
    self.message("Adjoint sensitivities without parameter or control quadratures")
    x = ca.MX.sym("x", 3)
    rates = ca.DM([1, 2, 3])
    x0 = ca.MX.sym("x0", 3)
    for backend in ["cvodes", "idas"]:
      if not ca.has_integrator(backend):
        continue
      for grid in [[1], [0.5, 1]]:
        I = ca.integrator("I", backend, {"x": x, "ode": -rates*x}, 0, grid,
          {"enable_forward": False, "abstol": 1e-11, "reltol": 1e-11})
        xf = I(x0=x0)["xf"]
        J = ca.Function("J", [x0], [ca.jacobian(xf[:2, -1], x0)])
        reference = ca.diag(ca.exp(-rates))[:2, :]
        for value in [[1, 2, 3], [3, 2, 1]]:
          self.checkarray(reference, J(value), digits=7)

  @requires_integrator('idas')
  def test_constraints_idas(self):
    x = ca.SX.sym("x")
    z = ca.SX.sym("z")
    p = ca.SX.sym("p")

    dae = { 'x':x, 'z':z, 'p':p, 'ode':z+p, 'alg':z*ca.cos(z)-x }

    opts = {"constraints": [1, 1]}
    I = ca.integrator("I","idas", dae, opts)
    sol = I(x0=0, p=0.15)
    # xf:0.259754>=0, zf:0.26948>=0

    opts = {"constraints": [-1, -1]}
    I = ca.integrator("I","idas", dae, opts)

    with self.assertInException("IDA_CONSTR_FAIL"):
      sol = I(x0=0, p=0.15)
      # xf:0.259754<=0, zf:0.26948<=0

  @requires_integrator('idas')
  def test_constraints_idas_sensitivities(self):
    self.message("IDAS constraints apply to primal states, not sensitivities")
    x = ca.MX.sym("x", 2)
    z = ca.MX.sym("z")
    p = ca.MX.sym("p", 2)
    I = ca.integrator("I", "idas",
      {"x": x, "z": z, "p": p, "ode": -p*x, "alg": z+ca.sum1(x)}, 0, 1,
      {"constraints": [2, 2, -2], "enable_reverse": False,
       "abstol": 1e-11, "reltol": 1e-11})
    result = I(x0=[1, 2], z0=-3, p=p)
    y = ca.vertcat(result["xf"], result["zf"])
    J = ca.Function("J", [p], [ca.jacobian(y, p)])
    reference = ca.DM([[-np.exp(-1), 0], [0, -2*np.exp(-2)],
                       [np.exp(-1), 2*np.exp(-2)]])
    self.checkarray(reference, J([1, 2]), digits=7)
    self.checkarray(reference, ca.Function.deserialize(J.serialize())([1, 2]), digits=7)

  @requires_integrator('idas')
  @requires_nlpsol('ipopt')
  def test_reduce_index(self):
    
    def problems():
      free  = 0 # default
      force = -1
      approx = 1
      blank = 0.0 # Will not be enforced but still used as initial guess

      # 2D Pendulum
      x = ca.SX.sym("x") # horizontal position of point
      y = ca.SX.sym("y") # vertical position of point

      dx = ca.SX.sym("dx")
      dy = ca.SX.sym("dy")

      u = ca.SX.sym("u") # helper states to bring second order system to first order
      v = ca.SX.sym("v")

      du = ca.SX.sym("du")
      dv = ca.SX.sym("dv")

      # Force
      T = ca.SX.sym("T")


      L = 1
      g = 9.81

      alg = ca.vertcat(dx-u,dy-v,du-T*x,dv-T*y+g,x**2+y**2-L**2)

      # Implicit differential states
      x_impl  = ca.vertcat(x,y,u,v)
      # Derivatives of implict differential states
      dx_impl = ca.vertcat(dx,dy,du,dv)

      z = T


      init_strength = {}
      init_strength["x_impl"] = ca.vertcat(free,force,free,free)
      init_strength["dx_impl"] = ca.vertcat(force,free,free,free)
      init_strength["z"] = ca.vertcat(free)
      init = {}
      init["x_impl"]  = ca.vertcat(-1.0,-0.5,blank,blank)  
      init["dx_impl"] = ca.vertcat(-0.1,blank,blank,blank)
      init["z"]  = ca.vertcat(blank)

      dae = {"x_impl": x_impl, "dx_impl": dx_impl, "z": z, "alg": alg}

      yield (dae,3,init_strength,init,5e-6)

      # 3D Pendulum


      x = ca.SX.sym("x")
      y = ca.SX.sym("y")
      z = ca.SX.sym("z")

      vx = ca.SX.sym("vx")
      vy = ca.SX.sym("vy")
      vz = ca.SX.sym("vz")

      dx = ca.SX.sym("dx")
      dy = ca.SX.sym("dy")
      dz = ca.SX.sym("dz")

      dvx = ca.SX.sym("dvx")
      dvy = ca.SX.sym("dvy")
      dvz = ca.SX.sym("dvz")

      lambd = ca.SX.sym("lambd")


      L = 1
      g = 9.81

      alg = ca.vertcat(dx-vx,dy-vy,dz-vz,dvx-lambd*x,dvy-lambd*y+g,dvz-lambd*z,x**2+y**2+z**2-L**2)

      x_impl = ca.vertcat(x,y,z,vx,vy,vz)
      dx_impl = ca.vertcat(dx,dy,dz,dvx,dvy,dvz)

      z = lambd

      init_strength = {}
      init_strength["x_impl"] = ca.vertcat(free,force,force,free,free,free)
      init_strength["dx_impl"] = ca.vertcat(force,force,free,free,free,free)
      init_strength["z"] = ca.vertcat(free)
      init = {}
      init["x_impl"]  = ca.vertcat(-1.0,-0.5,-0.5,blank,blank,blank)  
      init["dx_impl"] = ca.vertcat(-0.1,0.2,blank,blank,blank,blank)
      init["z"]  = ca.vertcat(blank)

      dae = {"x_impl": x_impl, "dx_impl": dx_impl, "z": z, "alg": alg}

      yield (dae,3,init_strength,init,1e-5)


      ## Flash system https://github.com/joaoleal/CppADCodeGen/blob/v2.4.3/test/cppad/cg/dae_index_reduction/model/flash.hpp 
      nEthanol = ca.SX.sym("nEthanol")
      nWater = ca.SX.sym("nWater")
      T = ca.SX.sym("T")
      yWater = ca.SX.sym("yWater")
      yEthanol = ca.SX.sym("yEthanol")
      FV = ca.SX.sym("FV")

      Q = 500
      F_feed = 10
      p = 1
      xFEthanol = 0.5
      T_feed = 50


      dnEthanol = ca.SX.sym("dnEthanol")
      dnWater = ca.SX.sym("dnWater")
      dT = ca.SX.sym("dT")

      alg = []

      FNEthanol = xFEthanol * F_feed
      n_lWater = nWater * 1000.
      n_lEthanol = nEthanol * 1000.
      m_lEthanol = n_lEthanol * 0.04606844
      m_lWater = n_lWater * 0.0180152833
      m = m_lEthanol + m_lWater
      wEthanol = m_lEthanol / m
      w_Water = 1 - wEthanol
      rho = 1 / (w_Water / 983.159471259596 + wEthanol / 743.278841365274)
      V = m / rho
      cWater = n_lWater / V
      dh = V / 0.502654824574367
      p_aux = p * 101325.
      dp = 101325. - p_aux
      v = ca.sqrt((9.80665 * dh + dp / rho) * 2.)
      F_VL = v * 0.0005
      FNLWater = cWater * F_VL
      n = n_lEthanol + n_lWater
      xWater = n_lWater / n
      F_NL = FNLWater / xWater
      FNLEthanol = F_NL - FNLWater
      FNVWater = yWater * FV
      FNVEthanol = FV - FNVWater
      D__nEthanol__Dt = FNEthanol - FNLEthanol - FNVEthanol
      alg.append(dnEthanol - (D__nEthanol__Dt * 0.001))

      FNFWater = F_feed - FNEthanol
      D__nWater__Dt = FNFWater - FNLWater - FNVWater
      alg.append(dnWater - (D__nWater__Dt * 0.001))

      FMFEthanol = FNEthanol * 0.04606844
      FMFWater = FNFWater * 0.0180152833
      Tfeed = T_feed + 273.15
      T_aux = T + 273.15
      dQ_F = (FMFEthanol * 2898.42878374706 + FMFWater * 4186.92536027523) * (Tfeed - T_aux)
      dH_mVapEthanol = 50430. * pow(1 - T_aux / 514., 0.4989) * ca.exp((0.4475 * T_aux) / 514.)
      T_r = T_aux / 647.
      dH_mVapWater = 52053. * pow(1 - T_r, 0.3199 + -0.212 * T_r + 0.25795 * T_r * T_r)
      dH_mVap = yEthanol * dH_mVapEthanol + yWater * dH_mVapWater
      dQ_vap = FV * dH_mVap
      Q_1 = Q * 1000.
      cp = wEthanol * 2898.42878374706 + w_Water * 4186.92536027523
      alg.append(dT - (dQ_F - dQ_vap + Q_1) / (m * cp))

      pVapWater = 100000. * pow(10., 4.6543 - 1435.264 / (T_aux - 64.848))
      KWater = pVapWater / p_aux
      alg.append(yWater - xWater * KWater)

      xEthanol = 1 - xWater
      pVapEthanol = ca.exp(74.475 + -7164.3 / T_aux + -7.327 * ca.log(T_aux) + 3.134e-06 * pow(T_aux, 2.))
      KEthanol = pVapEthanol / p_aux
      alg.append(yEthanol - xEthanol * KEthanol)
      
      alg.append(yWater + yEthanol - 1)


      alg = ca.vvcat(alg)

      x_impl = ca.vvcat([nEthanol,nWater,T])
      dx_impl = ca.vvcat([dnEthanol,dnWater,dT])

      z = ca.vertcat(yWater,yEthanol,FV)

      dae = {"x_impl": x_impl, "dx_impl": dx_impl, "z": z, "alg": alg}

      init_strength = {}
      init_strength["x_impl"] = ca.vertcat(force,force,approx)
      init_strength["dx_impl"] = ca.vertcat(free,free,free)
      init_strength["z"] = ca.vertcat(free,free,free)
      init = {}
      init["x_impl"]  = ca.vertcat(2.5,6.4,91)
      init["dx_impl"] = ca.vertcat(0.1,blank,blank)
      init["z"]  = ca.vertcat(0.5,0.5,0.5)

      dae = {"x_impl": x_impl, "dx_impl": dx_impl, "z": z, "alg": alg}

      yield (dae,2,init_strength,init,1e-8)

      # Exothermic Reactor from pantelides paper


      C  = ca.SX.sym("C")
      dC = ca.SX.sym("dC")
      T  = ca.SX.sym("T")
      dT = ca.SX.sym("dT")
      R  = ca.SX.sym("R")

      C0 = 1
      T0 = 1
      K1 = K2 = K3 = K4 = 1
      Tc = 1

      alg = ca.vertcat(K1*(C0-C)-R-dC,K1*(T0-T)+K2*R-K3*(T-Tc)-dT,R-K3*ca.exp(-K4/T)*C)

      x_impl = ca.vvcat([C,T])
      dx_impl = ca.vvcat([dC,dT])

      z = ca.vertcat(R)

      dae = {"x_impl": x_impl, "dx_impl": dx_impl, "z": z, "alg": alg}

      init_strength = {}
      init_strength["x_impl"] = ca.vertcat(free,force)
      init_strength["dx_impl"] = ca.vertcat(free,free)
      init_strength["z"] = ca.vertcat(free)
      init = {}
      init["x_impl"]  = ca.vertcat(blank,2)
      init["dx_impl"] = ca.vertcat(blank,blank)
      init["z"]  = ca.vertcat(blank)

      yield (dae,1,init_strength,init,1e-8)


      # S. Mattsson and G. Soederlind, Index reduction in differential-algebraic equations using dummy derivatives

      x1 = ca.SX.sym("x1")
      x2 = ca.SX.sym("x2")
      x3 = ca.SX.sym("x3")
      x4 = ca.SX.sym("x4")
      x5 = ca.SX.sym("x5")
      x6 = ca.SX.sym("x6")
      x7 = ca.SX.sym("x7")

      dx1 = ca.SX.sym("dx1")
      dx2 = ca.SX.sym("dx2")
      dx3 = ca.SX.sym("dx3")
      dx4 = ca.SX.sym("dx4")
      dx5 = ca.SX.sym("dx5")
      dx6 = ca.SX.sym("dx6")
      dx7 = ca.SX.sym("dx7")

      eqs = ca.vvcat([dx1+dx2,x2-5])
      u1 = 1
      u2 = 1
      u3 = 1
      u4 = 1

      x_impl = ca.vvcat([x1,x2,x3,x4,x5,x6,x7])
      dx_impl = ca.vvcat([dx1,dx2,dx3,dx4,dx5,dx6,dx7])

      eqs = ca.vcat([x5-dx1,x6-dx2,x7-dx3,x1+x2+u1,x1+x2+x3+u2,dx3+x1+x4+u3,2*dx5 + dx6 + dx7 + dx4 + u4])
      z = ca.SX(0,1)

      dae = {"x_impl": x_impl, "dx_impl": dx_impl, "z": z, "alg": eqs}

      # 2 for soares secchi
      yield (dae,3,None,None,None)

      # 

      """
      x1 = SX.sym("x1")
      x2 = SX.sym("x2")

      dx1 = SX.sym("dx1")
      dx2 = SX.sym("dx2")

      eqs = vvcat([dx1+dx2,x2-5])

      x_impl = vvcat([x1,x2])
      dx_impl = vvcat([dx1,dx2])

      z = SX(0,1)

      dae = {"x_impl": x_impl, "dx_impl": dx_impl, "z": z, "alg": eqs}
      (dae_reduced,constraints,var_map) = reduce_index(dae,gamma=-1)

      print(dae_reduced)
      print(constraints)


      x_impl0 = vertcat(-1,1)
      x_impl0_strength = vertcat(2,0) # strength 1: fix, 0: initial guess
      
      dx_impl0 = vertcat(0,0)
      dx_impl0_strength = vertcat(0,0) # strength 1: fix, 0: initial guess

      z0 = DM(0,1)
      z0_strength = DM(0,1) # strength 1: fix, 0: initial guess


      initial = {"x_impl": x_impl0,"dx_impl": dx_impl0,"z": z0}
      initial_strength = {"x_impl": x_impl0_strength,"dx_impl": dx_impl0_strength,"z": z0_strength}
      print("integrate")
      integrate(dae_reduced,constraints,var_map,initial,initial_strength)
      """

      # Uncontrollable DAE Pantelides

      x  = ca.SX.sym("x")
      dx = ca.SX.sym("dx")

      u1 = ca.SX.sym("u1")
      u2 = ca.SX.sym("u2")

      t = ca.SX.sym("t")

      y1 = t
      y2 = ca.sin(t)


      alg = ca.vertcat(x+u1+u2,x+dx+y1,x+y2)

      # Implicit differential states
      x_impl  = ca.vertcat(x)
      # Derivatives of implict differential states
      dx_impl = ca.vertcat(dx)

      z = ca.vertcat(u1,u2)

      dae = {"x_impl": x_impl, "dx_impl": dx_impl, "z": z, "alg": alg, "t": t}

      yield (dae,"Structural detection failed",None,None,None)


    print("loop")
    for [dae_orig,index,init_strength_orig,init_orig,max_rms] in problems():
      print("dae",dae_orig)
      for perturb_seed in range(5):
        if perturb_seed==0:
          dae = dae_orig
          init_strength = init_strength_orig
          init = init_orig
        else:
          np.random.seed(perturb_seed)

          dae = copy.copy(dae_orig)
          init_strength = copy.copy(init_strength_orig)
          init = copy.copy(init_orig)
          for k in ["x_impl","z","alg"]:
            s = ca.vertsplit(dae[k])
            if len(s)==0: continue
            perm = np.random.permutation(range(len(s)))
            dae[k] = ca.vcat([s[e] for e in perm])  # pyright: ignore[reportArgumentType]
            if init is not None and k in init: init[k] = init[k][perm]
            if init_strength is not None and k in init_strength: init_strength[k] = init_strength[k][perm]
            if k=="x_impl":
              s = ca.vertsplit(dae["d"+k])
              dae["d"+k] = ca.vcat([s[e] for e in perm])  # pyright: ignore[reportArgumentType]
              if init is not None: init["d"+k] = init["d"+k][perm]
              if init_strength is not None: init_strength["d"+k] = init_strength["d"+k][perm]


        grid = list(np.linspace(0,5,1000))
        if isinstance(index,str):
          with self.assertInException(index):
            ca.dae_reduce_index(dae, {"baumgarte_pole": -1, "algorithm": "pantelides", "max_iter": 20})
          continue
        (dae_reduced,meta) = ca.dae_reduce_index(dae, {"baumgarte_pole": -1, "algorithm": "pantelides"}) # soares_secchi
        assert meta["index"]==index
        if init_strength is not None:
          [dae_se, state_to_orig, phi] = ca.dae_map_semi_expl(dae, dae_reduced)
          intg = ca.integrator("intg","idas",dae_se,0, grid)
          init_gen = ca.dae_init_gen(dae, dae_reduced, "ipopt", init_strength, {"error_on_fail" : True})
          xz0 = init_gen(**init)  # pyright: ignore[reportCallIssue]
          sol = intg(**xz0)
          error = phi(x=sol["xf"],z=sol["zf"])["I"]
          #import pylab as plt
          #plt.plot(error.T)
          #plt.show()
      
          rms = np.sqrt(np.mean(error**2,axis=1))
          print("rms error", rms)
          self.assertTrue(np.all(rms<=max_rms))

  def test_IDACalcIC_u(self):
    x = ca.MX.sym("x")
    u = ca.MX.sym("u")
    z = ca.MX.sym("z")

    ivp = {}
    ivp["x"] = x
    ivp["z"] = z
    ivp["u"] = u
    ivp["ode"] = z
    ivp["alg"] = ca.sqrt(u)-z

    intg = ca.integrator("intg","idas",ivp,0,1)

    res = intg(x0=0,u=1)
    self.assertAlmostEqual(res["xf"],1,5,"test_IDACalcIC_u")

  @requires_integrator("cvodes")
  def test_expand(self):
    t=ca.MX.sym("t")
    q=ca.MX.sym("q")
    p=ca.MX.sym("p")

    dae = {'t':t, 'x':q, 'p':p, 'ode':q/p*t**2}
    integrator = ca.integrator("integrator", "cvodes", dae, 0, 2.3, {"expand":True})
    
    self.assertTrue(integrator.get_function('jacF').is_a("SXFunction"))
    
    t=ca.MX.sym("t")
    q=ca.MX.sym("q")
    p=ca.MX.sym("p")

    dae = {'t':t, 'x':q, 'p':p, 'ode':q/p*t**2}
    integrator = ca.integrator("integrator", "cvodes", dae, 0, 2.3, {"expand":False})
    
    self.assertTrue(integrator.get_function('jacF').is_a("MXFunction"))
    
  @requires_integrator("cvodes")
  def test_postpone_expand(self):
    x = ca.MX.sym("x")
    p = ca.MX.sym("p")
    
    J = ca.Function("jac_f", [x, ca.MX(1,1)], [-x], ['x', 'out_r'], ['jac_r_x'],{"always_inline":True})
    f = ca.Function('f', [x], [x**2], ['x'], ['r'], dict(custom_jacobian = J, jac_penalty = 0))
    
    dae = {'x':x, 'ode':f(x-2)}
    
    integrator = ca.integrator("integrator", "cvodes", dae, 0, 0.1, {"expand":False})
    
    integrator(x0=1)
    self.assertTrue(integrator.get_function('jacF').is_a("MXFunction"))
    self.checkarray(integrator.get_function('jacF')(x=1)["jac_ode_x"],1)
    
    dae = {'x':x, 'ode':f(x-2)}
    
    integrator = ca.integrator("integrator", "cvodes", dae, 0, 0.1, {"expand":True})
    integrator(x0=1)
    self.assertTrue(integrator.get_function('jacF').is_a("SXFunction"))
    self.checkarray(integrator.get_function('jacF')(x=1)["jac_ode_x"],-2)
    
    integrator = ca.integrator("integrator", "cvodes", dae, 0, 0.1, {"expand":False, "postpone_expand":True})
    
    integrator(x0=1)
    self.assertTrue(integrator.get_function('jacF').is_a("MXFunction"))
    self.checkarray(integrator.get_function('jacF')(x=1)["jac_ode_x"],1)
    
    dae = {'x':x, 'ode':f(x-2)}
    
    integrator = ca.integrator("integrator", "cvodes", dae, 0, 0.1, {"expand":True, "postpone_expand":True})
    integrator(x0=1)
    self.assertTrue(integrator.get_function('jacF').is_a("SXFunction"))
    self.checkarray(integrator.get_function('jacF')(x=1)["jac_ode_x"],1)

  @requires_integrator("idas")
  def test_idas_adjoint_zf_p(self):
    # issue #4063: direct dependence of zf on p and u in reverse mode
    x = ca.MX.sym('x'); z = ca.MX.sym('z'); p = ca.MX.sym('p'); u = ca.MX.sym('u')
    for alg in [z-p*x-u**2, z-p*x]:
      for ode, grid in [(-x, [1]), (-x+u, [1])]:
        dae = {'x':x,'z':z,'p':p,'u':u,'ode':ode,'alg':alg}
        F = ca.integrator('F','idas',dae,0,grid,{'abstol':1e-12,'reltol':1e-12})
        P = ca.MX.sym('P'); U = ca.MX.sym('U',1,len(grid))
        out = F(x0=2,z0=0,p=P,u=U)
        e = ca.vertcat(ca.vec(out['zf']), ca.vec(out['xf']))
        w = ca.vertcat(P,ca.vec(U))
        J_fwd = ca.Function('J_fwd',[w],[ca.jacobian(e,w,{"helper_options":{"ad_weight":0}})])
        J_rev = ca.Function('J_rev',[w],[ca.jacobian(e,w,{"helper_options":{"ad_weight":1}})])
        v = ca.DM([0.5,0.3,0.7][:1+len(grid)])
        self.checkarray(J_fwd(v), J_rev(v), digits=6)

  @requires_integrator("cvodes")
  def test_jit_fwd_factory(self):
    # issue #4061
    if not args.run_slow: return
    x = ca.MX.sym('x')
    dae = {'x': x, 'ode': -x}
    opts = {'abstol':1e-10,'reltol':1e-10,'jit':True,'compiler':'shell'}
    F = ca.integrator('F','cvodes',dae,0,[0.5,1],opts)
    self.assertTrue(F.get_function('jacF').is_a('External'))
    I_fwd = F.factory('I_fwd',['x0','fwd:x0'],['xf','fwd:xf'])
    [xf, fwd_xf] = I_fwd(2.0,1.0)
    self.checkarray(xf, 2*np.exp(-ca.DM([[0.5,1]])), digits=7)
    self.checkarray(fwd_xf, np.exp(-ca.DM([[0.5,1]])), digits=7)

  def test_issue3371(self):
    # Regression test for https://github.com/casadi/casadi/issues/3371
    # FixedStepIntegrator::advance_noevent / impulseB used the trailing nz_
    # entries of the augmented v / rv vectors, which lands inside a single
    # forward-seed block. The result was that fwd_zf (forward) and the
    # forward sensitivities of the adjoint outputs (forward-of-reverse)
    # were wrong for batched forward(N>1) calls.
    X = ca.MX.sym("x", 2)
    U = ca.MX.sym("u")
    Z = ca.MX.sym("z")
    dae = {"x": X, "p": U, "z": Z,
           "ode": ca.vertcat(Z * X[0] - X[1] + U, X[0]),
           "alg": Z - (1 - X[1]**2)}

    def cross_check(I, base_inputs, fwd_seeds_per_seed, msg):
      """Check batched forward(N) matches N independent forward(1) calls."""
      N = len(fwd_seeds_per_seed)
      I1 = I.forward(1)
      IN = I.forward(N)
      fwd_names = [n for n in IN.name_in() if n.startswith("fwd_")]
      # Build batched inputs by horzcat'ing per-seed entries
      batched_seeds = {}
      for n in fwd_names:
        cols = []
        for s in fwd_seeds_per_seed:
          v = s.get(n, ca.DM.zeros(I1.sparsity_in(n)))
          cols.append(v)
        batched_seeds[n] = ca.horzcat(*cols)
      batched_call = dict(base_inputs)
      batched_call.update(batched_seeds)
      batched = IN.call(batched_call)
      per_seed = {}
      for s in fwd_seeds_per_seed:
        single = {n: s.get(n, ca.DM.zeros(I1.sparsity_in(n))) for n in fwd_names}
        single_call = dict(base_inputs)
        single_call.update(single)
        out = I1.call(single_call)
        for k, v in out.items():
          per_seed.setdefault(k, []).append(v)
      for k, b in batched.items():
        merged = ca.horzcat(*per_seed[k])
        self.checkarray(b, merged, msg + " :: " + k)

    # --- Forward sensitivity: trips line 2074 (advance_noevent) ---
    I = ca.integrator("I", "collocation", dae)
    base = {"x0": ca.DM([0.5, 0.0]), "z0": ca.DM([1.0]), "p": ca.DM([0.1]),
            "u": ca.DM(0, 1)}
    base["out_xf"], base["out_zf"], base["out_qf"] = (
      I(x0=base["x0"], z0=base["z0"], p=base["p"])[k]
      for k in ("xf", "zf", "qf"))
    seeds = [
      {"fwd_x0": ca.DM([1.0, 0.0])},
      {"fwd_x0": ca.DM([0.0, 1.0])},
      {"fwd_p": ca.DM([1.0])},
    ]
    cross_check(I, base, seeds, "fwd-of-collocation")

    # Pin the smaller reproducer from issue #3371 (comment 1): with x0 = p = 0
    # the trajectory stays at x=0, z=1, so fwd_zf is zero for both seeds.
    base0 = {"x0": ca.DM.zeros(2), "z0": ca.DM([1.0]), "p": ca.DM([0.0]),
             "u": ca.DM(0, 1),
             "out_xf": ca.DM.zeros(2, 1), "out_zf": ca.DM([1.0]),
             "out_qf": ca.DM(0, 1)}
    I2 = I.forward(2)
    base0_call = dict(base0)
    base0_call.update({"fwd_x0": ca.DM([[0, 1], [0, 0]]),
                       "fwd_z0": ca.DM.zeros(I2.sparsity_in("fwd_z0")),
                       "fwd_p":  ca.DM([1, 0]).T,
                       "fwd_u":  ca.DM.zeros(I2.sparsity_in("fwd_u"))})
    out = I2.call(base0_call)
    self.checkarray(out["fwd_zf"], ca.DM.zeros(1, 2), "issue3371 fwd_zf")

    # --- Adjoint sensitivity: trips line 2251 (impulseB) ---
    # Forward-of-reverse: the augmented integrator has nfwd_>0 AND nadj_>0.
    Iadj = ca.integrator("I", "collocation", dae, {"nadj": 1})
    ref = Iadj(x0=ca.DM([0.5, 0.0]), z0=ca.DM([1.0]), p=ca.DM([0.1]),
               adj_xf=ca.DM([0.0, 0.0]), adj_zf=ca.DM([1.0]),
               adj_qf=ca.DM(0, 1))
    base = {"x0": ca.DM([0.5, 0.0]), "z0": ca.DM([1.0]), "p": ca.DM([0.1]),
            "u": ca.DM(0, 1),
            "adj_xf": ca.DM([0.0, 0.0]), "adj_zf": ca.DM([1.0]),
            "adj_qf": ca.DM(0, 1),
            "out_xf": ref["xf"], "out_zf": ref["zf"], "out_qf": ref["qf"],
            "out_adj_x0": ref["adj_x0"], "out_adj_z0": ref["adj_z0"],
            "out_adj_p": ref["adj_p"], "out_adj_u": ref["adj_u"]}
    seeds = [
      {"fwd_adj_zf": ca.DM([1.0])},
      {"fwd_x0": ca.DM([1.0, 0.0])},
      {"fwd_p": ca.DM([1.0])},
    ]
    cross_check(Iadj, base, seeds, "fwd-of-reverse-of-collocation")

  def check_jacobian(self, out, v, v0, digits):
    # Forward and adjoint Jacobians of out(v) at v0 agree
    J_fwd = ca.Function("J", [v], [ca.jtimes(out, v, ca.MX.eye(v.numel()))])(v0)
    J_adj = ca.Function("J", [v], [ca.jtimes(out, v, ca.MX.eye(out.numel()), True).T])(v0)
    self.checkarray(J_fwd, J_adj, "fwd vs adj", digits=digits)

  def check_hessian(self, obj, v, v0, digits):
    # Hessian of obj(v) at v0: forward and adjoint mode over the forward and adjoint gradients
    n = v.numel()
    g_fwd = ca.jtimes(obj, v, ca.MX.eye(n)).T
    g_adj = ca.jtimes(obj, v, ca.MX.ones(1, 1), True)
    H = [ca.Function("H", [v], [h])(v0) for h in [ca.jtimes(g_fwd, v, ca.MX.eye(n)),
         ca.jtimes(g_adj, v, ca.MX.eye(n)), ca.jtimes(g_fwd, v, ca.MX.eye(n), True).T,
         ca.jtimes(g_adj, v, ca.MX.eye(n), True).T]]
    for h in H[1:]: self.checkarray(H[0], h, "Hessian", digits=digits)
    return H[0]

  def test_reverse_grid_directions(self):
    self.message("Several adjoint directions with several output times")
    x = ca.SX.sym("x", 2)
    p = ca.SX.sym("p")
    dae = {"x": x, "p": p, "ode": ca.vertcat(x[1], -p*x[0]**3), "quad": x[0]**2}
    for plugin, opts in [("cvodes", {"reltol": 1e-10, "abstol": 1e-10}),
                         ("idas", {"reltol": 1e-10, "abstol": 1e-10}),
                         ("rosenbrock", {"step_size": 0.125}),
                         ("erk", {"method": "rk45", "reltol": 1e-8, "abstol": 1e-8})]:
      if not ca.has_integrator(plugin): continue
      F = ca.integrator("F", plugin, dae, 0, [0.5, 1.0], opts)
      inputs = [ca.DM([0.3, -0.2]), ca.DM(), ca.DM(0.7), ca.DM(1, 0)]
      inputs += [ca.DM(F.sparsity_in(i)) for i in range(4, F.n_in())]
      outputs = list(F(*inputs))
      seeds = [[ca.DM(F.sparsity_out(i), numpy.random.default_rng(10*d+i).random(F.nnz_out(i)))
                for i in range(F.n_out())] for d in range(2)]
      adj2 = F.reverse(2)(*(inputs + outputs + [ca.horzcat(seeds[0][i], seeds[1][i])
                                               for i in range(F.n_out())]))
      for d in range(2):
        adj1 = F.reverse(1)(*(inputs + outputs + seeds[d]))
        for i in range(F.n_in()):
          if F.nnz_in(i) == 0: continue
          self.checkarray(adj1[i], ca.horzsplit_n(adj2[i], 2)[d], plugin, digits=8)

  @memory_heavy()
  def test_events_sparsity(self):
    self.message("zero-crossing events: dependencies through the transition")
    g = 9.81
    x = ca.SX.sym("x", 2)
    p = ca.SX.sym("p")
    # Bouncing ball: the restitution p enters the transition only, which also couples h and v
    dae = {"x": x, "p": p, "ode": ca.vertcat(x[1], -g), "zero": x[0]}
    grid = list(numpy.linspace(0.1, 4, 40))
    for plugin, opts in [("rk", {"number_of_finite_elements": 100}),
                         ("cvodes", {"abstol": 1e-12, "reltol": 1e-12})]:
      if not ca.has_integrator(plugin): continue
      F = ca.integrator("F", plugin, dae, 0, grid, dict(opts, transition=bounce_transition(),
        event_tol=1e-12, max_event_iter=20))
      v = ca.MX.sym("v", 3)
      xf = F(x0=v[:2], p=v[2])["xf"][:, -1]
      G = ca.Function("G", [v], [xf, ca.jacobian(xf, v)])
      self.assertTrue(G.sparsity_out(1).is_dense())
      v0 = ca.DM([2, 0.3, 0.8])
      self.checkarray(G(v0)[1], fd_jacobian(G, v0), plugin, digits=5)

  @memory_heavy()
  @requires_integrator("idas")
  def test_events_idas(self):
    self.message("idas: zero-crossing events, also of a DAE")
    g = 9.81
    x = ca.SX.sym("x", 2)
    z = ca.SX.sym("z")
    p = ca.SX.sym("p")
    # A single bounce from height 2 at rest, at sqrt(2*2/g), restitution p
    grid = list(numpy.linspace(0.1, 1.5, 15))
    tb = numpy.sqrt(4/g)
    ub = 0.8*g*tb
    h = [2 - g*tt**2/2 if tt < tb else ub*(tt - tb) - g*(tt - tb)**2/2 for tt in grid]
    v = [-g*tt if tt < tb else ub - g*(tt - tb) for tt in grid]
    ode = {"x": x, "p": p, "ode": ca.vertcat(x[1], -g), "zero": x[0]}
    dae = {"x": x, "z": z, "p": p, "ode": ca.vertcat(z, -g), "alg": z - x[1], "zero": x[0]}
    for problem, transition in [(ode, bounce_transition()), (dae, bounce_transition(1))]:
      F = ca.integrator("F", "idas", problem, 0, grid, {"transition": transition,
        "abstol": 1e-12, "reltol": 1e-12, "event_tol": 1e-12, "max_event_iter": 20})
      r = F(x0=[2, 0], z0=0, p=0.8)
      self.checkarray(r["xf"], ca.DM([h, v]), "xf", digits=7)

  def test_reverse_over_adjoint_directions(self):
    self.message("Derivatives of adjoint sensitivities with several directions (#3293)")
    x = ca.SX.sym("x", 3)
    p = ca.SX.sym("p", 2)
    dae = {"x": x, "p": p, "ode": ca.vertcat(-p[0]*x[0], -p[1]*x[1]*x[0], -p[1]*x[2]**2),
           "quad": x[0]*x[2]}
    rng = numpy.random.default_rng(1)
    rand = lambda sp: ca.DM(sp, rng.random(sp.nnz()))
    for plugin, opts, digits in [("cvodes", {"reltol": 1e-12, "abstol": 1e-12}, 7),
                                 ("idas", {"reltol": 1e-12, "abstol": 1e-12}, 7),
                                 ("rosenbrock", {"step_size": 1./6}, 12),
                                 ("rosenbrock", {"reltol": 1e-6, "abstol": 1e-6}, 11),
                                 ("erk", {"step_size": 0.2}, 12),
                                 ("erk", {"method": "rk45", "reltol": 1e-6, "abstol": 1e-6}, 11)]:
      if not ca.has_integrator(plugin): continue
      for grid in [[1.0], [0.5, 1.0]]:
        F = ca.integrator("F", plugin, dae, 0, grid, opts)
        # Adjoint sensitivities with two directions, and their derivatives
        B = F.reverse(2)
        inputs = [rand(B.sparsity_in(i)) for i in range(B.n_in())]
        outputs = list(B(*inputs))
        v = [rand(B.sparsity_in(i)) for i in range(B.n_in())]
        w = [[rand(B.sparsity_out(i)) for i in range(B.n_out())] for d in range(2)]
        fwd = B.forward(1)(*(inputs + outputs + v))
        adj2 = B.reverse(2)(*(inputs + outputs + [ca.horzcat(w[0][i], w[1][i])
                                                 for i in range(B.n_out())]))
        for d in range(2):
          adj = B.reverse(1)(*(inputs + outputs + w[d]))
          self.checkarray(sum(float(ca.dot(w[d][i], fwd[i])) for i in range(B.n_out())),
                          sum(float(ca.dot(adj[i], v[i])) for i in range(B.n_in())),
                          plugin + " duality", digits=digits)
          for i in range(B.n_in()):
            if B.nnz_in(i) > 0:
              self.checkarray(adj[i], ca.horzsplit_n(adj2[i], 2)[d], plugin, digits=digits)
        # Jacobian, relying on the sparsity pattern, vs dense forward directional derivatives
        vin = ca.MX.sym("v", sum(B.nnz_in(i) for i in range(B.n_in())))
        args = ca.vertsplit(vin, [int(e) for e in numpy.cumsum([0] + [B.nnz_in(i) for i in range(B.n_in())])])
        res = B(*[ca.MX(B.sparsity_in(i), args[i]) for i in range(B.n_in())])
        vout = ca.vertcat(*[ca.vec(r.nz[:]) for r in res])
        J = ca.Function("J", [vin], [ca.jacobian(vout, vin),
                                     ca.jtimes(vout, vin, ca.MX.eye(vin.numel()))])
        J_sp, J_dense = J(numpy.concatenate([e.nonzeros() for e in inputs]))
        self.checkarray(J_dense, J_sp, plugin + " sparsity", digits=digits)

  def test_runtime_integrator_quadrature_error_control(self):
    # Constant states cannot reveal error in a time-varying quadrature.
    t = ca.SX.sym("t")
    x = ca.SX.sym("x")
    dae = {"t": t, "x": x, "ode": 0*x, "quad": ca.exp(10*t)}
    exact = numpy.expm1(10)/10
    opts = {"abstol": 1e-9, "reltol": 1e-9, "step0": 1.0}
    for plugin in ["erk", "rosenbrock"]:
      if not ca.has_integrator(plugin):
        continue
      default = ca.integrator("default", plugin, dae, 0, 1, opts)
      enabled = ca.integrator("enabled", plugin, dae, 0, 1,
                              dict(opts, quad_err_con=True))
      disabled = ca.integrator("disabled", plugin, dae, 0, 1,
                               dict(opts, quad_err_con=False))
      q_default = float(default(x0=0)["qf"])
      q_enabled = float(enabled(x0=0)["qf"])
      q_disabled = float(disabled(x0=0)["qf"])
      self.assertEqual(q_default, q_enabled)
      self.assertLess(abs(q_default-exact), 1e-4)
      self.assertGreater(abs(q_disabled-exact), 1e-2)
      self.assertGreater(default.stats()["n_steps"], disabled.stats()["n_steps"])

  @memory_heavy()
  def test_rosenbrock_order(self):
    self.message("rosenbrock: order 5 on fixed meshes, ODE and index-1 DAE")
    t = ca.SX.sym("t")
    x = ca.SX.sym("x", 2)
    z = ca.SX.sym("z")
    dae = {"t": t, "x": x, "z": z, "ode": ca.vertcat(-z + 0.5*x[1], -x[0]*z + ca.sin(t)),
           "alg": z**3 + z - x[0]**2*(1 + 0.5*ca.cos(t)), "quad": z*x[1]}
    ref = ca.integrator("ref", "idas", dae, 0, 2, {"reltol": 1e-14, "abstol": 1e-14})(x0=[1, 0.5], z0=0.6)
    err = []
    for N in [8, 16, 32]:
      F = ca.integrator("F", "rosenbrock", dae, 0, 2, {"step_size": 2./N})
      r = F(x0=[1, 0.5], z0=0.6)
      err.append(max(float(ca.norm_inf(r[k] - ref[k])) for k in ["xf", "zf", "qf"]))
    for e0, e1 in zip(err[:-1], err[1:]):
      self.assertGreater(numpy.log2(e0/e1), 4.6)

  @memory_heavy()
  def test_rosenbrock_ind(self):
    self.message("rosenbrock: derivatives are exact derivatives of the scheme")
    t = ca.SX.sym("t")
    x = ca.SX.sym("x", 2)
    z = ca.SX.sym("z")
    p = ca.SX.sym("p", 2)
    u = ca.SX.sym("u")
    dae = {"t": t, "x": x, "z": z, "p": p, "u": u,
           "ode": ca.vertcat(-z + p[0]*x[1] + u, -x[0]*z + p[1]*ca.sin(t)),
           "alg": z**3 + z - x[0]**2*(1 + 0.5*ca.cos(t)), "quad": z*x[1] + p[1]*x[0]**2}
    # Including an output at t0, with a seed on the consistent algebraic state there
    for opts, grid in [(o, g) for o in [{"step_size": 0.1}, {"reltol": 1e-6, "abstol": 1e-6}]
                       for g in [[0.5, 1, 2], [0, 1, 2]]]:
      F = ca.integrator("F", "rosenbrock", dae, 0, grid, opts)
      v = ca.MX.sym("v", 7)
      r = F(x0=v[:2], z0=0.5, p=v[2:4], u=v[4:].T)
      out = ca.vertcat(ca.vec(r["xf"]), ca.vec(r["zf"]), ca.vec(r["qf"]))
      obj = ca.sumsqr(out)
      v0 = [1, 0.5, 0.5, 0.3, 0.1, -0.2, 0.3]
      # Forward and adjoint sensitivities agree to machine precision, to second order
      self.check_jacobian(out, v, v0, 12)
      H = self.check_hessian(obj, v, v0, 12)
      self.checkarray(H, H.T, "Hessian symmetry", digits=12)

  @memory_heavy()
  def test_rosenbrock_dense_output(self):
    self.message("rosenbrock: outputs inside steps")
    t = ca.SX.sym("t")
    x = ca.SX.sym("x", 2)
    z = ca.SX.sym("z")
    p = ca.SX.sym("p", 2)
    u = ca.SX.sym("u")
    dae = {"t": t, "x": x, "z": z, "p": p, "u": u,
           "ode": ca.vertcat(-z + p[0]*x[1] + u, -x[0]*z + p[1]*ca.sin(t)),
           "alg": z**3 + z - x[0]**2*(1 + 0.5*ca.cos(t)) - u, "quad": z*x[1] + p[1]*x[0]**2}
    ode = {"t": t, "x": x, "p": p, "u": u, "ode": ca.vertcat(x[1] + p[0]*u, -p[1]*x[0] + ca.sin(t)),
           "quad": x[0]*x[1]**2}
    for sys, ref_plugin in [(ode, "cvodes"), (dae, "idas")]:
      if not ca.has_integrator(ref_plugin): continue
      args = {"x0": [1, 0.5], "p": [0.5, 0.3]}
      # Consistent z0 for the reference: z^3 + z = 1.5 x0^2 - u
      if "z" in sys: args["z0"] = float(min(numpy.roots([1, 0, 1, -1.3]), key=lambda r: abs(r.imag)).real)
      for opts in [{"step_size": 0.04}, {"reltol": 1e-6, "abstol": 1e-6}]:
        # Within a control group, the number of steps does not depend on the number of outputs,
        # accuracy inside steps
        n_steps = []
        for grid in [[2.], list(numpy.linspace(0.002, 2, 1000))]:
          F = ca.integrator("F", "rosenbrock", sys, 0, grid,
                            dict(opts, control_group_size=[len(grid)]))
          r = F(u=0.2, **args)
          n_steps.append(F.stats()["n_steps"])
          ref = ca.integrator("ref", ref_plugin, sys, 0, grid,
            {"reltol": 1e-10, "abstol": 1e-10})(u=0.2, **args)
          for k in ["xf", "zf", "qf"]:
            self.checkarray(ref[k], r[k], k, digits=5)
        self.assertTrue(n_steps[1] <= n_steps[0] + 2)
        # Derivatives are exact, with constant controls and with a step in the controls
        grid = [0.13, 0.31, 0.37, 0.52, 0.9, 1.15, 1.2, 1.71, 2.0]
        for u0 in [[0.1] * 9, [0.1] * 4 + [0.3] * 5]:
          F = ca.integrator("F", "rosenbrock", sys, 0, grid, opts)
          v = ca.MX.sym("v", 13)
          r = F(x0=v[:2], z0=0.5, p=v[2:4], u=v[4:].T) if "z" in sys else F(x0=v[:2], p=v[2:4],
            u=v[4:].T)
          out = ca.vertcat(ca.vec(r["xf"]), ca.vec(r["zf"]), ca.vec(r["qf"]))
          v0 = [1, 0.5, 0.5, 0.3] + u0
          self.check_jacobian(out, v, v0, 12)
          self.check_hessian(ca.sumsqr(out), v, v0, 10)

  @memory_heavy()
  def test_rosenbrock_codegen(self):
    self.message("rosenbrock: code generation, including derivatives")
    x = ca.SX.sym("x", 2)
    z = ca.SX.sym("z")
    p = ca.SX.sym("p")
    dae = {"x": x, "z": z, "p": p, "ode": ca.vertcat(x[1], -p*z), "alg": z**3 + z - x[0],
           "quad": x[0]**2}
    for opts, grid in [({"reltol": 1e-8, "abstol": 1e-8, "max_tape_steps": 1000}, [0.5, 1]),
                       ({"step_size": 0.05}, [0.5, 1]),
                       ({"step_size": 0.03}, list(numpy.linspace(0.07, 1, 15))),
                       ({"reltol": 1e-6, "abstol": 1e-6, "max_tape_steps": 1000},
                        list(numpy.linspace(0.07, 1, 15)))]:
      F = ca.integrator("F", "rosenbrock", dae, 0, grid, opts)
      v = ca.MX.sym("v", 3)
      r = F(x0=v[:2], p=v[2])
      obj = ca.sumsqr(r["xf"]) + ca.sumsqr(r["qf"])
      outs = [r["xf"], r["zf"], r["qf"], ca.gradient(obj, v)]
      # Forward over forward: steps replayed from a nominal integration
      if len(grid) > 2: outs.append(ca.jtimes(ca.jtimes(obj, v, ca.MX.eye(3)).T, v, ca.MX.eye(3)))
      G = ca.Function("G", [v], outs)
      self.check_codegen(G, inputs=[ca.DM([0.4, -0.3, 2])], digits=12)

  @memory_heavy()
  def test_erk(self):
    self.message("erk: order 4, dense output, steps independent of the outputs")
    t = ca.SX.sym("t")
    x = ca.SX.sym("x", 2)
    p = ca.SX.sym("p")
    u = ca.SX.sym("u")
    dae = {"t": t, "x": x, "p": p, "u": u, "ode": ca.vertcat(x[1], -p*x[0] + u*ca.cos(t)),
           "quad": x[0]**2}
    # Harmonic oscillator: error at the end and inside a step
    w = 1.3
    exact = lambda tt: ca.DM([numpy.cos(w*tt), -w*numpy.sin(w*tt)])
    err = []
    for h in [0.1, 0.05, 0.025]:
      r = ca.integrator("F", "erk", dae, 0, [0.37, 2.0], {"method": "rk4", "step_size": h})(x0=[1, 0],
        p=w**2)
      err.append([float(ca.norm_inf(r["xf"][:, 1] - exact(2.0))),
                  float(ca.norm_inf(r["xf"][:, 0] - exact(0.37)))])
    for e0, e1 in zip(err[:-1], err[1:]):
      self.assertGreater(numpy.log2(e0[0]/e1[0]), 3.8)
      self.assertGreater(numpy.log2(e0[1]/e1[1]), 3.4)
    # Same steps as rk when the output times are on the mesh
    args = {"x0": [1, 0.2], "p": 2, "u": ca.DM([[0.3, -0.1]])}
    r = ca.integrator("F", "erk", dae, 0, [0.5, 1.0], {"method": "rk4", "step_size": 0.1})(**args)
    ref = ca.integrator("F", "rk", dae, 0, [0.5, 1.0], {"number_of_finite_elements": 10})(**args)
    for k in ["xf", "qf"]: self.checkarray(ref[k], r[k], k, digits=14)
    # Within a control group, the number of steps does not depend on the number of outputs
    for n in [10, 1000]:
      F = ca.integrator("F", "erk", dae, 0, list(numpy.linspace(0, 10, n)),
                        {"step_size": 0.05, "control_group_size": [n]})
      F(x0=[1, 0], p=w**2)
      self.assertEqual(F.stats()["n_steps"], 200)

  @memory_heavy()
  def test_erk_derivatives(self):
    self.message("erk: derivatives are exact derivatives of the scheme")
    t = ca.SX.sym("t")
    x = ca.SX.sym("x", 2)
    p = ca.SX.sym("p")
    u = ca.SX.sym("u")
    dae = {"t": t, "x": x, "p": p, "u": u, "ode": ca.vertcat(x[1], -p*x[0]**3 + u*ca.cos(t)),
           "quad": x[0]**2 + u*x[1]}
    # Outputs inside steps, at t0, and controls that change or not
    grid = [0, 0.13, 0.31, 0.37, 1.0, 1.15, 2.0]
    for u0 in [[0.1] * 7, [0.1] * 3 + [0.3] * 4, list(numpy.linspace(0.1, 0.4, 7))]:
      F = ca.integrator("F", "erk", dae, 0, grid, {"step_size": 0.1})
      v = ca.MX.sym("v", 10)
      r = F(x0=v[:2], p=v[2], u=v[3:].T)
      out = ca.vertcat(ca.vec(r["xf"]), ca.vec(r["qf"]))
      v0 = [1, 0.5, 1.7] + u0
      self.check_jacobian(out, v, v0, 13)
      obj = ca.sumsqr(out)
      self.check_hessian(obj, v, v0, 12)
      # Finite differences, also for controls of equal value
      g = ca.Function("g", [v], [ca.gradient(obj, v)])(v0)
      self.checkarray(g, fd_jacobian(ca.Function("f", [v], [obj]), v0).T, "FD", digits=7)

  @memory_heavy()
  def test_erk_codegen(self):
    self.message("erk: code generation, including derivatives")
    x = ca.SX.sym("x", 2)
    p = ca.SX.sym("p")
    u = ca.SX.sym("u")
    dae = {"x": x, "p": p, "u": u, "ode": ca.vertcat(x[1], -p*ca.sin(x[0]) + u),
           "quad": x[0]**2}
    for grid, opts in [([0.5, 1], {"step_size": 0.04}),
                       (list(numpy.linspace(0.07, 1, 15)), {"step_size": 0.04}),
                       (list(numpy.linspace(0.07, 1, 15)), {"method": "rk4", "step_size": 0.1}),
                       (list(numpy.linspace(0.07, 1, 15)),
                        {"method": "rk45", "reltol": 1e-6, "abstol": 1e-6, "max_tape_steps": 200})]:
      F = ca.integrator("F", "erk", dae, 0, grid, opts)
      v = ca.MX.sym("v", 3 + len(grid))
      r = F(x0=v[:2], p=v[2], u=v[3:].T)
      obj = ca.sumsqr(r["xf"]) + ca.sumsqr(r["qf"])
      # Forward over forward in a few directions: all of them make a large file to compile
      S = ca.DM.eye(v.numel())[:, [0, 2, 3]]
      outs = [r["xf"], r["qf"], ca.gradient(obj, v), ca.jtimes(ca.jtimes(obj, v, S).T, v, S)]
      G = ca.Function("G", [v], outs)
      v0 = [0.4, -0.3, 2] + [0.1 * (i % 3) for i in range(len(grid))]
      self.check_codegen(G, inputs=[ca.DM(v0)], digits=12)

  @memory_heavy()
  def test_control_groups(self):
    self.message("erk, rosenbrock: one control per control group, derivatives for each control")
    t = ca.SX.sym("t")
    x = ca.SX.sym("x", 2)
    z = ca.SX.sym("z")
    p = ca.SX.sym("p")
    u = ca.SX.sym("u")
    ode = {"t": t, "x": x, "p": p, "u": u, "ode": ca.vertcat(x[1], -p*x[0]**3 + u*ca.cos(t)),
           "quad": x[0]**2 + u*x[1]}
    dae = {"t": t, "x": x, "z": z, "p": p, "u": u, "ode": ca.vertcat(x[1], -p*x[0] + u*z),
           "alg": z**3 + z - u - 0.5*x[0], "quad": x[0]**2 + u*z}
    # Output times on the mesh of fixed steps
    grid = [0.2, 0.3, 0.5, 1.0, 1.2, 2.0]
    for plugin, opts in [("erk", {"method": "rk4", "step_size": 0.1}),
                         ("rosenbrock", {"step_size": 0.1}),
                         ("erk", {"reltol": 1e-10, "abstol": 1e-10, "max_tape_steps": 500}),
                         ("rosenbrock", {"reltol": 1e-10, "abstol": 1e-10, "max_tape_steps": 500})]:
      d = dae if plugin == "rosenbrock" else ode
      # Repeated controls and their control groups agree, up to the restarts between outputs
      args = {"x0": [1, 0.5], "p": 1.7}
      r = ca.integrator("F", plugin, d, 0, grid, opts)(u=ca.DM([[0.1, 0.1, 0.3, 0.2, 0.2, 0.2]]),
                                                      **args)
      F = ca.integrator("F", plugin, d, 0, grid, dict(opts, control_group_size=[2, 1, 3]))
      r_group = F(u=ca.DM([[0.1, 0.3, 0.2]]), **args)
      for k in ["xf", "zf", "qf"]: self.checkarray(r[k], r_group[k], k, digits=6)
      # Controls of equal value: each has its own sensitivities
      for cgs in [None, [2, 1, 3]]:
        F = ca.integrator("F", plugin, d, 0, grid, dict(opts, control_group_size=cgs or [1] * 6))
        n = 3 + F.size2_in(3)
        v = ca.MX.sym("v", n)
        r = F(x0=v[:2], p=v[2], u=v[3:].T)
        out = ca.vertcat(ca.vec(r["xf"]), ca.vec(r["zf"]), ca.vec(r["qf"]))
        G = ca.Function("G", [v], [out, ca.jtimes(out, v, ca.MX.eye(n)),
                                   ca.jtimes(out, v, ca.MX.eye(out.numel()), True).T])
        v0 = numpy.array([1, 0.5, 1.7] + [0.2] * (n - 3))
        _, J_fwd, J_adj = G(v0)
        J_fd = fd_jacobian(G, v0)
        self.checkarray(J_fd, J_fwd, "fwd vs FD", digits=7)
        self.checkarray(J_fd, J_adj, "adj vs FD", digits=7)
        if cgs: self.check_codegen(G, inputs=[ca.DM(v0)], digits=12)
        # States depend on the controls of their group and earlier ones, forward and reverse
        start = numpy.cumsum([0] + (cgs or [1] * 6))
        for w in [0, 1]:
          J = ca.Function("J", [v], [out], {"ad_weight_sp": w}).jac_sparsity(0, 0)
          dep = numpy.array(ca.DM(J, 1))[:12, 3:] != 0
          for k in range(6):
            for j in range(n - 3):
              self.assertEqual(bool(dep[2*k:2*k+2, j].any()), bool(start[j] <= k))
    with self.assertInException("must sum to the number of output times"):
      ca.integrator("F", "erk", ode, 0, grid, {"control_group_size": [2, 1]})

  @memory_heavy()
  def test_erk_rk45(self):
    self.message("erk: pair of order 5(4), fixed and adaptive steps")
    t = ca.SX.sym("t")
    x = ca.SX.sym("x", 2)
    p = ca.SX.sym("p")
    dae = {"t": t, "x": x, "p": p, "ode": ca.vertcat(x[1], -p*x[0]), "quad": x[0]**2}
    w = 1.3
    exact = lambda tt: numpy.array([numpy.cos(w*tt), -w*numpy.sin(w*tt)])
    # Fixed steps: order 5
    err = []
    for h in [0.1, 0.05, 0.025]:
      r = ca.integrator("F", "erk", dae, 0, 2, {"method": "rk45", "step_size": h})(x0=[1, 0], p=w**2)
      err.append(numpy.max(numpy.abs(numpy.array(r["xf"]).ravel() - exact(2.0))))
    for e0, e1 in zip(err[:-1], err[1:]):
      self.assertGreater(numpy.log2(e0/e1), 4.7)
    # Adaptive steps: error follows the tolerance; the last stage is the first of the next step
    grid = list(numpy.linspace(0.1, 10, 50))
    for tol in [1e-6, 1e-9]:
      F = ca.integrator("F", "erk", dae, 0, grid, {"method": "rk45", "reltol": tol, "abstol": tol})
      r = F(x0=[1, 0], p=w**2)
      e = numpy.max(numpy.abs(numpy.array(r["xf"]) - numpy.array([exact(tt) for tt in grid]).T))
      self.assertLess(e, 10*tol)
      self.assertGreater(e, 0.01*tol)
      st = F.stats()
      self.assertEqual(st["n_rhs"], 1 + 6*(st["n_steps"] + st["n_reject"]))
    # Rejected steps: a fast transient with a large initial step
    dae2 = {"x": x, "ode": ca.vertcat(-50*(x[0] - ca.cos(x[1])), 1)}
    F = ca.integrator("F", "erk", dae2, 0, 1, {"method": "rk45", "step0": 0.5})
    r = F(x0=[0, 0])
    ref = ca.integrator("F", "cvodes", dae2, 0, 1, {"reltol": 1e-12, "abstol": 1e-12})(x0=[0, 0])
    self.checkarray(ref["xf"], r["xf"], "rejections", digits=5)
    self.assertGreater(F.stats()["n_reject"], 0)
    with self.assertInException("unknown method"):
      ca.integrator("F", "erk", dae, 0, 1, {"method": "rk5"})

  @memory_heavy()
  def test_erk_rk45_derivatives(self):
    self.message("erk: adaptive steps, derivatives of the scheme on the nominal mesh")
    t = ca.SX.sym("t")
    x = ca.SX.sym("x", 2)
    p = ca.SX.sym("p")
    u = ca.SX.sym("u")
    dae = {"t": t, "x": x, "p": p, "u": u, "ode": ca.vertcat(x[1], -p*x[0]**3 + u*ca.cos(t)),
           "quad": x[0]**2 + u*x[1]}
    grid = [0, 0.13, 0.31, 0.37, 1.0, 1.15, 2.0]
    for u0 in [[0.1] * 7, [0.1] * 3 + [0.3] * 4]:
      F = ca.integrator("F", "erk", dae, 0, grid, {"method": "rk45", "reltol": 1e-5,
                                                   "abstol": 1e-5})
      v = ca.MX.sym("v", 10)
      r = F(x0=v[:2], p=v[2], u=v[3:].T)
      out = ca.vertcat(ca.vec(r["xf"]), ca.vec(r["qf"]))
      v0 = [1, 0.5, 1.7] + u0
      self.check_jacobian(out, v, v0, 13)
      # Forward over forward replays the nominal steps
      self.check_hessian(ca.sumsqr(out), v, v0, 11)

  @memory_heavy()
  def test_events_bouncing_ball(self):
    self.message("zero-crossing events: values, forward sensitivities, code generation")
    g = 9.81
    x = ca.SX.sym("x", 2)
    p = ca.SX.sym("p")
    dae = {"x": x, "p": p, "ode": ca.vertcat(x[1], -g), "quad": x[0], "zero": x[0]}
    transition = bounce_transition()
    def exact(h0, p0, tt):
      # Height, velocity and integral of the height
      tb, ub, q = 0., None, 0.
      h, v = h0, 0.
      while True:
        tau = numpy.sqrt(2*h/g) if ub is None else 2*ub/g
        if tb + tau > tt: break
        q += h*tau - g*tau**3/6 if ub is None else ub*tau**2/2 - g*tau**3/6
        v = (-g*tau if ub is None else ub - g*tau)
        tb, ub, h = tb + tau, -p0*v, 0.
      tau = tt - tb
      u0 = 0. if ub is None else ub
      return [h + u0*tau - g*tau**2/2, u0 - g*tau], q + h*tau + u0*tau**2/2 - g*tau**3/6
    grid = list(numpy.linspace(0.1, 4, 40))
    h0, p0 = 2., 0.8
    xf_ref = numpy.array([exact(h0, p0, tt)[0] for tt in grid]).T
    qf_ref = numpy.array([exact(h0, p0, tt)[1] for tt in grid])
    for plugin, opts, features in event_integrators:
      if not features["events"]: continue
      F = ca.integrator("F", plugin, dae, 0, grid, dict(opts, transition=transition,
        event_tol=1e-12))
      r = F(x0=[h0, 0], p=p0)
      self.checkarray(r["xf"], ca.DM(xf_ref), plugin, digits=features["digits"])
      self.checkarray(r["qf"], ca.DM(qf_ref).T, plugin, digits=features["digits"])
      if "n_events" in F.stats(): self.assertEqual(F.stats()["n_events"], 5)
      # Forward sensitivities through the events, against finite differences
      v = ca.MX.sym("v", 3)
      rv = F(x0=v[:2], p=v[2])
      out = ca.vertcat(rv["xf"][:, -1], rv["qf"][-1])
      G = ca.Function("G", [v], [out, ca.jacobian(out, v)])
      v0 = ca.DM([h0, 0.3, p0])
      self.checkarray(G(v0)[1], fd_jacobian(G, v0), plugin, digits=min(6, features["digits"]))
      if features["codegen"]: self.check_codegen(G, inputs=[v0], digits=10)

  @memory_heavy()
  def test_events_multiple(self):
    self.message("zero-crossing events: several functions, transition by event index")
    g, ceiling = 9.81, 1.5
    x = ca.SX.sym("x", 2)
    p = ca.SX.sym("p")
    index = ca.SX.sym("index")
    t = ca.SX.sym("t")
    # Ball between a floor (event 0, restitution p) and a ceiling (event 1, elastic)
    dae = {"x": x, "p": p, "ode": ca.vertcat(x[1], -g), "zero": ca.vertcat(x[0], ceiling - x[0])}
    transition = ca.Function("transition", [index, t, x, ca.SX(0, 1), p, ca.SX(0, 1)],
      [ca.vertcat(x[0], ca.if_else(index == 0, -p, -1)*x[1]), ca.SX(0, 1)],
      ["index", "t", "x", "z", "p", "u"], ["post_x", "post_z"])
    def exact(v0, tt):
      # Piecewise parabolas: the first crossing of the floor or the ceiling ends each one
      t0, h, v = 0., float(v0[0]), float(v0[1])
      while True:
        taus = []
        for level, sign in [(0., -1), (ceiling, 1)]:
          # h + v tau - g tau^2 / 2 = level, crossing in the direction sign
          disc = v*v - 2*g*(level - h)
          if disc < 0: continue
          for tau in [(v - numpy.sqrt(disc))/g, (v + numpy.sqrt(disc))/g]:
            if tau > 1e-12 and sign*(v - g*tau) > 0: taus.append((tau, level))
        tau, level = min(taus)
        if t0 + tau > tt: break
        v = (v - g*tau)*(-float(v0[2]) if level == 0 else -1)
        t0, h = t0 + tau, level
      tau = tt - t0
      return [h + v*tau - g*tau**2/2, v - g*tau]
    grid = list(numpy.linspace(0.1, 3, 30))
    v0 = ca.DM([1, 4, 0.8])
    xf_ref = ca.DM(numpy.array([exact(v0, tt) for tt in grid]).T)
    for plugin, opts, features in event_integrators:
      if not features["events"]: continue
      F = ca.integrator("F", plugin, dae, 0, grid, dict(opts, transition=transition,
        event_tol=1e-12))
      v = ca.MX.sym("v", 3)
      xf = F(x0=v[:2], p=v[2])["xf"]
      G = ca.Function("G", [v], [ca.vec(xf), ca.jacobian(ca.vec(xf), v)])
      self.checkarray(G(v0)[0], ca.vec(xf_ref), plugin, digits=features["digits"])
      self.checkarray(G(v0)[1], fd_jacobian(G, v0), plugin, digits=min(6, features["digits"]))
      if features["codegen"]: self.check_codegen(G, inputs=[v0], digits=10)

  @memory_heavy()
  def test_events_dae(self):
    self.message("zero-crossing events of a DAE, algebraic states consistent after the transition")
    g = 9.81
    x = ca.SX.sym("x", 2)
    z = ca.SX.sym("z")
    p = ca.SX.sym("p")
    dae = {"x": x, "z": z, "p": p, "ode": ca.vertcat(z, -g), "alg": z - x[1], "zero": x[0]}
    # The transition leaves z inconsistent
    transition = bounce_transition(1)
    # Several bounces, outputs between an event and the end of the step that contains it
    grid = list(numpy.linspace(0.1, 3, 30))
    ref = ca.integrator("ref", "erk", {"x": x, "p": p, "ode": ca.vertcat(x[1], -g), "zero": x[0]},
      0, grid, {"transition": bounce_transition(), "event_tol": 1e-12, "abstol": 1e-11,
      "reltol": 1e-11})
    v0 = ca.DM([2, 0.3, 0.8])
    xf_ref = ref(x0=v0[:2], p=v0[2])["xf"]
    for plugin, opts, features in event_integrators:
      if not features["events_dae"]: continue
      F = ca.integrator("F", plugin, dae, 0, grid, dict(opts, transition=transition,
        event_tol=1e-12))
      v = ca.MX.sym("v", 3)
      r = F(x0=v[:2], z0=v[1], p=v[2])
      out = ca.vertcat(ca.vec(r["xf"]), ca.vec(r["zf"]))
      G = ca.Function("G", [v], [out, ca.jacobian(out, v)])
      self.checkarray(G(v0)[0], ca.vertcat(ca.vec(xf_ref), xf_ref[1, :].T), plugin,
        digits=features["digits"])
      # Sensitivities through the events, with the rates after the transition
      self.checkarray(G(v0)[1], fd_jacobian(G, v0), plugin, digits=min(6, features["digits"]))
      if features["codegen"]: self.check_codegen(G, inputs=[v0], digits=10)

  @memory_heavy()
  def test_events_breaking_spring(self):
    self.message("zero-crossing events with controls and parameters")
    x = ca.SX.sym("x", 3)
    p = ca.SX.sym("p", 3)
    u = ca.SX.sym("u")
    index = ca.SX.sym("index")
    t = ca.SX.sym("t")
    # States v, position, broken; parameters mass, stiffness, damping; control disturbance
    f = ca.if_else(x[2] < 0.5, -p[1]*x[1] + u, 0)
    dae = {"x": x, "p": p, "u": u, "ode": ca.vertcat((f - p[2]*x[0])/p[0], x[0], 0),
           "zero": 2 - x[1]}
    transition = ca.Function("transition", [index, t, x, ca.SX(0, 1), p, u],
      [ca.vertcat(x[0], x[1], 1), ca.SX(0, 1)], ["index", "t", "x", "z", "p", "u"],
      ["post_x", "post_z"])
    grid = list(numpy.linspace(0.04, 4, 100))
    d = numpy.zeros(100)
    d[20:30] = -5
    d[40:50] = 15
    x0, p0 = [-5, -1, 0], ca.DM([1, 2, 0.1])
    v = ca.MX.sym("v", 3)
    def outputs(F):
      xf = F(x0=x0, p=v, u=d)["xf"]
      return ca.Function("G", [v], [xf, ca.jacobian(xf, v)])
    ref = outputs(ca.integrator("ref", "erk", dae, 0, grid, {"transition": transition,
      "event_tol": 1e-12, "abstol": 1e-12, "reltol": 1e-12}))(p0)
    for plugin, opts, features in event_integrators:
      if not features["events"]: continue
      F = ca.integrator("F", plugin, dae, 0, grid, dict(opts, transition=transition,
        event_tol=1e-12))
      G = outputs(F)
      r = G(p0)
      self.checkarray(r[0], ref[0], plugin, digits=features["digits"] - 1)
      self.checkarray(r[1], ref[1], plugin, digits=features["digits"] - 2)
      F(x0=x0, p=p0, u=d)
      if "n_events" in F.stats(): self.assertEqual(F.stats()["n_events"], 1)
      if features["codegen"]: self.check_codegen(G, inputs=[p0], digits=12)
    # Every event counts towards max_events
    F = ca.integrator("F", "erk", dae, 0, grid, {"transition": transition, "max_events": 0})
    with self.assertInException("max_events"):
      F(x0=x0, p=p0, u=d)

  @memory_heavy()
  def test_events_adjoint(self):
    self.message("zero-crossing events: adjoint and forward-over-adjoint sensitivities")
    g = 9.81
    x = ca.SX.sym("x", 2)
    z = ca.SX.sym("z")
    p = ca.SX.sym("p")
    # Bouncing ball with restitution p, as an ODE and as a DAE with z the velocity
    ode = {"x": x, "p": p, "ode": ca.vertcat(x[1], -g), "quad": x[0], "zero": x[0]}
    dae = {"x": x, "z": z, "p": p, "ode": ca.vertcat(z, -g), "alg": z - x[1], "quad": x[0],
           "zero": x[0]}
    # Outputs before, at and after the events, also inside steps cut short by them
    grid = list(numpy.linspace(0.1, 3, 30))
    wts = ca.DM(numpy.random.default_rng(0).uniform(0.5, 1.5, (2, 30)))
    for plugin, opts, features in event_integrators:
      if not features["events_adjoint"]: continue
      for problem, transition in [(ode, bounce_transition())] + ([(dae, bounce_transition(1))]
          if features["events_dae"] else []):
        F = ca.integrator("F", plugin, problem, 0, grid, dict(opts, transition=transition,
          event_tol=1e-12))
        v = ca.MX.sym("v", 3)
        r = F(x0=v[:2], z0=v[1], p=v[2])
        obj = ca.sum1(ca.sum2(wts * r["xf"])) + ca.sumsqr(r["qf"])
        G = ca.Function("G", [v], [ca.jtimes(obj, v, ca.MX.ones(1, 1), True),
          ca.jtimes(obj, v, ca.MX.eye(3)).T, ca.hessian(obj, v)[0]])
        v0 = ca.DM([2, 0.3, 0.8])
        # Adjoint as forward sensitivities
        grad_adj, grad_fwd, H = G(v0)
        self.checkarray(grad_adj, grad_fwd, plugin, digits=9)
        # Forward-over-adjoint: central differences of the gradient, Richardson extrapolated
        fd = lambda h: fd_jacobian(G, v0, h)
        self.checkarray(H, (4*fd(1e-4) - fd(2e-4))/3, plugin, digits=6)
        if features["codegen"]: self.check_codegen(G, inputs=[v0], digits=9)

  @memory_heavy()
  def test_erk_tape(self):
    self.message("erk: adaptive steps and adjoints, bounded tape")
    x = ca.SX.sym("x", 2)
    dae = {"x": x, "ode": ca.vertcat(x[1], -x[0])}
    for n_tape in [0, 3, 1000]:
      F = ca.integrator("F", "erk", dae, 0, 10, {"method": "rk45", "max_tape_steps": n_tape})
      x0 = ca.MX.sym("x0", 2)
      G = ca.Function("G", [x0], [ca.jtimes(ca.sumsqr(F(x0=x0)["xf"]), x0, ca.MX.ones(1, 1), True)])
      if n_tape == 3:
        # More steps than the tape holds
        with self.assertInException("increase max_tape_steps"):
          G([1, 0.5])
      else:
        self.checkarray(G([1, 0.5]), ca.DM([2, 1]), "adjoint", digits=5)
      if n_tape == 0:
        # The tape grows in the virtual machine only
        with self.assertInException("need max_tape_steps for code generation"):
          ca.CodeGenerator("tape_cg.c").add(G)

if __name__ == '__main__':
    unittest.main()
