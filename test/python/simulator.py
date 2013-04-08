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

scipy_available = True
try:
	import scipy.special
	from scipy.linalg import expm
except:
	scipy_available = False

class Simulatortests(casadiTestCase):


  def setUp(self):
    # Reference solution is q0 e^((t^3-t0^3)/(3 p))
    t=ssym("t")
    q=ssym("q")
    p=ssym("p")
    f=SXFunction(daeIn(t=t,x=q,p=p),daeOut(ode=q/p*t**2))
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
    
  def test_sim_inputs(self):
    self.message("Simulator inputs")
    num = self.num
    tc = DMatrix(n.linspace(0,num['tend'],100))
    
    t=ssym("t")
    q=ssym("q")
    p=ssym("p")
    
    out = SXFunction(daeIn(t=t, x=q, p=p),[q,t,p])
    out.init()
        
    f=SXFunction(daeIn(t=t, x=q, p=p),daeOut(ode=q/p*t**2))
    f.init()
    integrator = CVodesIntegrator(f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    integrator.setOption("fsens_err_con", True)
    #integrator.setOption("verbose",True)
    integrator.setOption("t0",0)
    integrator.setOption("tf",2.3)
    integrator.init()
    sim = Simulator(integrator,out,tc)
    sim.init()
    sim.input(0).set([num['q0']])
    sim.input(1).set([num['p']])
    sim.evaluate()

    self.checkarray(sim.output(),num['q0']*exp(tc**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.output(1),tc,"Evaluation output mismatch")
    self.checkarray(sim.output(2),DMatrix.ones(tc.shape)*num['p'],"Evaluation output mismatch")
    
    f=SXFunction(daeIn(t=t, x=q, p=p),daeOut(ode=q/p*t**2))
    f.init()
    integrator = CVodesIntegrator(f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    integrator.setOption("fsens_err_con", True)
    #integrator.setOption("verbose",True)
    integrator.setOption("t0",0)
    integrator.setOption("tf",2.3)
    integrator.init()
    sim = Simulator(integrator,out,tc)
    sim.init()
    sim.input(0).set([num['q0']])
    sim.input(1).set([num['p']])
    sim.evaluate()

    self.checkarray(sim.output(),num['q0']*exp(tc**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.output(1),tc,"Evaluation output mismatch")
    self.checkarray(sim.output(2),DMatrix.ones(tc.shape)*num['p'],"Evaluation output mismatch")
    
    out = SXFunction(daeIn(t=t, x=q),[q,t])
    out.init()
    
    f=SXFunction(daeIn(t=t, x=q),daeOut(ode=q/num['p']*t**2))
    f.init()
    integrator = CVodesIntegrator(f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    integrator.setOption("fsens_err_con", True)
    #integrator.setOption("verbose",True)
    integrator.setOption("t0",0)
    integrator.setOption("tf",2.3)
    integrator.init()
    sim = Simulator(integrator,out,tc)
    sim.init()
    sim.input(0).set([num['q0']])
    sim.evaluate()

    self.checkarray(sim.output(),num['q0']*exp(tc**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.output(1),tc,"Evaluation output mismatch")
    
    f=SXFunction(daeIn(x=q),daeOut(ode=-q))
    f.init()
    integrator = CVodesIntegrator(f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    integrator.setOption("fsens_err_con", True)
    #integrator.setOption("verbose",True)
    integrator.setOption("t0",0)
    integrator.setOption("tf",2.3)
    integrator.init()
    sim = Simulator(integrator,out,tc)
    sim.init()
    sim.input(0).set([num['q0']])
    sim.evaluate()

    self.checkarray(sim.output(),num['q0']*exp(-tc),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.output(1),tc,"Evaluation output mismatch",digits=9)
    self.assertTrue(sim.output().sparsity()==sim.output(1).sparsity())
    
  def test_sim_outputs(self):
    self.message("Simulator: outputs")
    num=self.num
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    t = n.linspace(0,num['tend'],100)
    sim = Simulator(self.integrator,t)
    sim.init()
    sim.input(0).set([num['q0']])
    sim.input(1).set([num['p']])
    sim.evaluate()

    self.checkarray(sim.output(),DMatrix(q0*exp(t**3/(3*p))),"Evaluation output mismatch",digits=9)

    tv = SX("t")
    out = SXFunction(daeIn(t=tv),[tv])
    
    out.init()
    sim = Simulator(self.integrator,out,t)
    sim.init()
    sim.input(0).set([num['q0']])
    sim.input(1).set([num['p']])
    sim.evaluate()

    self.checkarray(sim.output(),t,"Evaluation output mismatch")

    pv = SX("p")
    out = SXFunction(daeIn(p=pv),[pv])
    
    out.init()
    sim = Simulator(self.integrator,out,t)
    sim.init()
    sim.input(0).set([num['q0']])
    sim.input(1).set([num['p']])
    sim.evaluate()

    
    self.checkarray(sim.output(),DMatrix.ones(sim.output().shape)*p,"Evaluation output mismatch")

    #yv = SX("y")
    #out = SXFunction(daeIn(),[yv])
    
    #out.init()
    #sim = Simulator(self.integrator,out,t)
    #sim.init()
    #sim.input(0).set([num['q0']])
    #sim.input(1).set([num['p']])
    #sim.evaluate()

    #self.checkarray(sim.output(),DMatrix.zeros(sim.output().shape),"INTEGRATOR_XPF unsupported",digits=9)
    #self.checkarray(sim.output(),DMatrix((q0*t**2*exp(t**3/(3*p)))/p),"Evaluation output mismatch",digits=9)

  def test_controlsim_inputs(self):
    self.message("ControlSimulator: inputs")
    num=self.num
    tc = 0.01*DMatrix([0,8,16,24,32])
    t  = ssym("t")
    q  = ssym("q")
    p  = ssym("p")
    
    t0 = ssym("t0")
    tf_= ssym("tf")
    
    qm = ssym("qm")
    
    out = SXFunction(controldaeIn(t=t, x=q, p=p, t0=t0, tf=tf_, x_major=qm),[q,t,p,t0,tf_,qm])
    out.init()
    
    f=SXFunction(controldaeIn(t=t, x=q, p=p, x_major=qm),[q/p*t**2])
    f.init()
    sim = ControlSimulator(f,out,tc)
    sim.setOption('nf',4)
    sim.setOption('integrator',CVodesIntegrator)
    sim.setOption('integrator_options', {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True})
    sim.init()
    sim.input("x0").set([num['q0']])
    sim.input("p").set([num['p']])
    self.assertTrue(sim.input("u").empty())
    sim.evaluate()
    
    tf = DMatrix(sim.getMinorT())


    self.checkarray(sim.output(),num['q0']*exp(tf**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.output(1),tf,"Evaluation output mismatch")
    self.checkarray(sim.output(2),DMatrix.ones(tf.shape)*num['p'],"Evaluation output mismatch")
    self.checkarray(sim.output(3),DMatrix([0,0,0,0,0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.24]),"Evaluation output mismatch")
    self.checkarray(sim.output(4),DMatrix([0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.32,0.32,0.32,0.32,0.32]),"Evaluation output mismatch")
    self.checkarray(sim.output(5),num['q0']*exp(sim.output(3)**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    
    
    #CasadiOptions.setCatchErrorsPython(False)
    f=SXFunction(controldaeIn(t=t, x=q, p=p, x_major=qm),[q/p*t**2])
    f.init()
    sim = ControlSimulator(f,out,tc)
    sim.setOption('nf',4)
    sim.setOption('integrator',CVodesIntegrator)
    sim.setOption('integrator_options', {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True})
    sim.init()
    sim.input("x0").set([num['q0']])
    sim.input("p").set([num['p']])
    sim.evaluate()
    
    tf = DMatrix(sim.getMinorT())

    self.checkarray(sim.output(),num['q0']*exp(tf**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.output(1),tf,"Evaluation output mismatch")
    self.checkarray(sim.output(2),DMatrix.ones(tf.shape)*num['p'],"Evaluation output mismatch")
    self.checkarray(sim.output(3),DMatrix([0,0,0,0,0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.24]),"Evaluation output mismatch")
    self.checkarray(sim.output(4),DMatrix([0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.32,0.32,0.32,0.32,0.32]),"Evaluation output mismatch")
    self.checkarray(sim.output(5),num['q0']*exp(sim.output(3)**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    
    out = SXFunction(controldaeIn(t=t, x=q, t0=t0, tf=tf_, x_major=qm),[q,t,t0,tf_,qm])
    out.init()
    
    f=SXFunction(controldaeIn(t=t, x=q, x_major=qm),[q/num['p']*t**2])
    f.init()
    sim = ControlSimulator(f,out,tc)
    sim.setOption('nf',4)
    sim.setOption('integrator',CVodesIntegrator)
    sim.setOption('integrator_options', {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True})
    sim.init()
    sim.input("x0").set([num['q0']])
    sim.evaluate()

    self.checkarray(sim.output(),num['q0']*exp(tf**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.output(1),tf,"Evaluation output mismatch")
    self.checkarray(sim.output(2),DMatrix([0,0,0,0,0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.24]),"Evaluation output mismatch")
    self.checkarray(sim.output(3),DMatrix([0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.32,0.32,0.32,0.32,0.32]),"Evaluation output mismatch")
    self.checkarray(sim.output(4),num['q0']*exp(sim.output(2)**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    
    f=SXFunction(controldaeIn(x=q, x_major=qm),[-q])
    f.init()
    sim = ControlSimulator(f,out,tc)
    sim.setOption('nf',4)
    sim.setOption('integrator',CVodesIntegrator)
    sim.setOption('integrator_options', {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True})
    sim.init()
    sim.input("x0").set([num['q0']])
    sim.evaluate()


    self.checkarray(sim.output(),num['q0']*exp(-tf),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.output(1),tf,"Evaluation output mismatch")
    self.checkarray(sim.output(2),DMatrix([0,0,0,0,0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.24]),"Evaluation output mismatch")
    self.checkarray(sim.output(3),DMatrix([0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.32,0.32,0.32,0.32,0.32]),"Evaluation output mismatch")
    self.checkarray(sim.output(4),num['q0']*exp(-sim.output(2)),"Evaluation output mismatch",digits=9)

  def test_controlsim_interpolation(self):
    q  = ssym("q")
    u  = ssym("u")
    ui = ssym("ui")


    
    tc = 0.01*DMatrix([0,8,16,24,32])
    
    U = DMatrix([0,0.1,0,0.2])
    
    q0=2.3

    out = SXFunction(controldaeIn(x=q, u=u, u_interp=ui),[q,u,ui])
    out.init()
    
    f=SXFunction(controldaeIn(x=q, u=u),[-q])
    f.init()
    sim = ControlSimulator(f,out,tc)
    sim.setOption('nf',4)
    sim.setOption('integrator',CVodesIntegrator)
    sim.setOption('integrator_options', {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True})
    self.assertRaises(Exception,lambda : sim.init())
    
    out = SXFunction(controldaeIn(x=q, u=u),[q,u])
    out.init()
    
    f=SXFunction(controldaeIn(x=q),[-q])
    f.init()
    sim = ControlSimulator(f,out,tc)
    sim.setOption('nf',4)
    sim.setOption('integrator',CVodesIntegrator)
    sim.setOption('integrator_options', {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True})
    self.assertRaises(Exception,lambda : sim.init())


    out = SXFunction(controldaeIn(x=q, u=u, u_interp=ui),[q,u,ui])
    out.init()
    
    f=SXFunction(controldaeIn(x=q, u=u, u_interp=ui),[-q])
    f.init()
    sim = ControlSimulator(f,out,tc)
    sim.setOption('nf',4)
    sim.setOption('integrator',CVodesIntegrator)
    sim.setOption('integrator_options', {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True})
    sim.init()
    sim.input("x0").set([2.3])
    sim.input("u").set(U)
    sim.evaluate()

    tf = DMatrix(sim.getMinorT())
    
    self.checkarray(sim.output(),q0*exp(-tf),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.output(1),DMatrix([0,0,0,0,0.1,0.1,0.1,0.1,0,0,0,0,0.2,0.2,0.2,0.2,0.2]),"Evaluation output mismatch")
    self.checkarray(sim.output(2),DMatrix([0,0.025,0.05,0.075,0.1,0.075,0.05,0.025,0,0.05,0.1,0.15,0.2,0.2,0.2,0.2,0.2]),"Evaluation output mismatch")

    out = SXFunction(controldaeIn(x=q, u=u, u_interp=ui),[q,u,ui])
    out.init()
    
    f=SXFunction(controldaeIn(x=q,u=u, u_interp=ui),[-q])
    f.init()
    sim = ControlSimulator(f,out,tc)
    sim.setOption('control_endpoint',True)
    sim.setOption('nf',4)
    sim.setOption('integrator',CVodesIntegrator)
    sim.setOption('integrator_options', {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True})
    sim.init()
    sim.input("x0").set([2.3])
    sim.input("u").set(vertcat([U,0]))
    sim.evaluate()

    tf = DMatrix(sim.getMinorT())
    
    self.checkarray(sim.output(),q0*exp(-tf),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.output(1),DMatrix([0,0,0,0,0.1,0.1,0.1,0.1,0,0,0,0,0.2,0.2,0.2,0.2,0.2]),"Evaluation output mismatch")
    self.checkarray(sim.output(2),DMatrix([0,0.025,0.05,0.075,0.1,0.075,0.05,0.025,0,0.05,0.1,0.15,0.2,0.15,0.1,0.05,0]),"Evaluation output mismatch")

  def test_controlsim_outputs(self):
    self.message("CVodes integration: outputs")
    num=self.num
    tc = 0.01*DMatrix([0,8,16,24,32])
    t=ssym("t")
    q=ssym("q")
    p=ssym("p")
    f=SXFunction(controldaeIn(t=t, x=q, p=p),[q/p*t**2])
    f.init()
    sim = ControlSimulator(f,tc)
    sim.setOption('nf',4)
    sim.setOption('integrator',CVodesIntegrator)
    sim.setOption('integrator_options', {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True})
    sim.init()
    sim.input(0).set([num['q0']])
    sim.input(1).set([num['p']])
    sim.evaluate()
    
    tf = DMatrix(sim.getMinorT())

    self.checkarray(sim.output(),num['q0']*exp(tf**3/(3*num['p'])),"Evaluation output mismatch",digits=9)


  def test_simulator_time_offset(self):
    self.message("CVodes integration: simulator time offset")
    num=self.num
    t = n.linspace(0.7,num['tend'],100)
    sim = Simulator(self.integrator,t)
    sim.init()
    sim.input(0).set([num['q0']])
    sim.input(1).set([num['p']])
    sim.evaluate()

    tend=num['tend']
    q0=num['q0']
    p=num['p']

    self.assertAlmostEqual(sim.output()[-1],q0*exp((tend**3-0.7**3)/(3*p)),9,"Evaluation output mismatch")
    
    
  def test_simulator_sensitivities(self):
    self.message("Forward sensitivities")
    t = SX("t")

    x  = SX("x") 
    v  = SX("v") 

    u = SX("u") 

    b = SX("b")
    b_ = 0.1
    k = SX("k")
    k_ = 1
    

    rhs = vertcat([v, ( -  b*v - k*x)])
    f=SXFunction(daeIn(t=t, x=vertcat([x,v]), p=vertcat([b,k])),daeOut(ode=rhs))
    f.init()
    
    cf=SXFunction(controldaeIn(t=t, x=vertcat([x,v]), p=vertcat([b,k])),[rhs])
    cf.init()

    x0 = SX("x0")
    dx0 = SX("dx0")


    X0 = DMatrix([1,0])
    p_ = [b_,k_]

    # algebraic solution
    sole = exp(-(b*t)/2)*((sin((sqrt(4*k-b**2)*t)/2)*(2*(dx0+x0*b)-x0*b))/sqrt(4*k-b**2)+x0*cos((sqrt(4*k-b**2)*t)/2))
    sol = SXFunction([t,vertcat([x0,dx0]),vertcat([b,k])],[vertcat([sole,jacobian(sole,t)])])
    sol.init()
    sol.input(0).set(50)
    sol.input(1).set(X0)
    sol.input(2).set(p_)


    for Integrator in [CVodesIntegrator, IdasIntegrator]:
      integrator = Integrator(f)
      #integrator.setOption("verbose",True)
      #integrator.setOption("monitor",["integrate"])
      integrator.setOption("fsens_abstol",1e-12)
      integrator.setOption("fsens_reltol",1e-12)
      integrator.setOption("fsens_err_con", True)
      integrator.setOption("tf",50)
      integrator.init()

      integrator.input("x0").set(X0)
      integrator.input("p").set(p_)
      integrator.fwdSeed("x0").set([1,0])
      integrator.fwdSeed("p").set([0,0])
      integrator.evaluate(1,0)
      
      fwdSens_int = DMatrix(integrator.fwdSens("xf"))

      ts = linspace(0,50,100)

      sim=Simulator(integrator,ts)
      sim.init()
      sim.input("x0").set(X0)
      sim.input("p").set(p_)
      sim.fwdSeed("x0").set([1,0])
      sim.fwdSeed("p").set([0,0])
      sim.evaluate(1,0)

      fwdSens_sim = DMatrix(sim.fwdSens("xf")[-1,:])
      
      csim = ControlSimulator(cf,ts)
      csim.setOption("integrator",Integrator)
      csim.setOption("integrator_options",{"fsens_abstol": 1e-12, "fsens_reltol": 1e-12, "fsens_err_con": True})
      csim.init()
      csim.input("x0").set(X0)
      csim.input("p").set(p_)
      csim.fwdSeed("x0").set([1,0])
      csim.fwdSeed("p").set([0,0])
      csim.evaluate(1,0)
      fwdSens_csim = DMatrix(csim.fwdSens("xf")[-1,:])
      
      sol.fwdSeed(1).set([1,0])
      sol.fwdSeed(2).set([0,0])
      sol.evaluate(1,0)
      fwdSens_exact = sol.fwdSens()
      
      #digits = 6
      #if (Integrator is IdasIntegrator):
      digits = 2 # Joel: No reason to treat Idas differently from CVodes

      self.assertAlmostEqual(fwdSens_int[0],fwdSens_exact[0],digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_int[1],fwdSens_exact[1],digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_sim[0],fwdSens_exact[0],digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_sim[1],fwdSens_exact[1],digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_csim[0],fwdSens_exact[0],digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_csim[1],fwdSens_exact[1],digits,"Forward sensitivity")
      
      integrator.fwdSeed("x0").set([0,0])
      integrator.fwdSeed("p").set([1,0])
      integrator.evaluate(1,0)
      
      fwdSens_int = DMatrix(integrator.fwdSens("xf"))
      
      sim.fwdSeed("x0").set([0,0])
      sim.fwdSeed("p").set([1,0])
      sim.evaluate(1,0)

      fwdSens_sim = DMatrix(sim.fwdSens("xf")[-1,:])

      csim.fwdSeed("x0").set([0,0])
      csim.fwdSeed("p").set([1,0])
      csim.evaluate(1,0)
      fwdSens_csim = DMatrix(sim.fwdSens("xf")[-1,:])
      
      sol.fwdSeed(1).set([0,0])
      sol.fwdSeed(2).set([1,0])
      sol.evaluate(1,0)
      fwdSens_exact = sol.fwdSens()
      
      #digits = 6
      #if (Integrator is IdasIntegrator):
      digits = 2 # Joel: No reason to treat Idas differently from CVodes

      self.assertAlmostEqual(fwdSens_int[0],fwdSens_exact[0], digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_int[1],fwdSens_exact[1], digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_sim[0],fwdSens_exact[0], digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_sim[1],fwdSens_exact[1], digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_csim[0],fwdSens_exact[0], digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_csim[1],fwdSens_exact[1], digits,"Forward sensitivity")

  def test_simulator_sensitivities_adj(self):
    # This test is currently disabled, awaiting support for the feature
    return
    
    self.message("Adjoint sensitivities")
    t = SX("t")

    x  = SX("x") 
    v  = SX("v") 

    u = SX("u") 

    b = SX("b")
    b_ = 0.1
    k = SX("k")
    k_ = 1
    

    rhs = vertcat([v, ( -  b*v - k*x)])
    f=SXFunction(daeIn(t=t, x=[x,v], p=[b,k]),daeOut(ode=rhs))
    f.init()
    
    cf=SXFunction(controldaeIn(t=t, x=[x,v], p=[b,k]),[rhs])
    cf.init()

    x0 = SX("x0")
    dx0 = SX("dx0")


    X0 = DMatrix([1,0])
    p_ = [b_,k_]
    
    
    Te = 50

    # algebraic solution
    sole = exp(-(b*t)/2)*((sin((sqrt(4*k-b**2)*t)/2)*(2*(dx0+x0*b)-x0*b))/sqrt(4*k-b**2)+x0*cos((sqrt(4*k-b**2)*t)/2))
    sol = SXFunction([[t],[x0,dx0],[b,k]],[vertcat([sole,jacobian(sole,t)])])
    sol.init()
    sol.input(0).set(Te)
    sol.input(1).set(X0)
    sol.input(2).set(p_)

    Integrator = CVodesIntegrator
  
    integrator = Integrator(f)
    #integrator.setOption("verbose",True)
    #integrator.setOption("monitor",["integrate"])
    integrator.setOption("abstolB",1e-12)
    integrator.setOption("reltolB",1e-12)
    integrator.setOption("fsens_err_con", True)
    integrator.setOption("tf",Te)
    integrator.init()

    integrator.input("x0").set(X0)
    integrator.input("p").set(p_)
    integrator.adjSeed("xf").set([1,0])
    integrator.evaluate(0,1)
    
    adjSens_X0_int = DMatrix(integrator.adjSens("x0"))
    adjSens_P_int = DMatrix(integrator.adjSens("p"))

    N = 100
    ts = linspace(0,Te,N)

    sim=Simulator(integrator,ts)
    sim.init()
    sim.input("x0").set(X0)
    sim.input("p").set(p_)
    sim.adjSeed(0).set([0,0]*(N-1) + [1,0])
    sim.evaluate(0,1)

    adjSens_X0_sim = DMatrix(sim.adjSens("x0"))
    adjSens_P_sim = DMatrix(sim.adjSens("p"))

    
    csim = ControlSimulator(cf,ts)
    csim.setOption("integrator",Integrator)
    csim.setOption("integrator_options",{"abstolB": 1e-12, "reltolB": 1e-12, "fsens_err_con": True, "steps_per_checkpoint":1000})
    csim.init()
    csim.input("x0").set(X0)
    csim.input("p").set(p_)
    csim.adjSeed(0).set([0,0]*(N-1) + [1,0])
    csim.evaluate(0,1)
    adjSens_X0_csim = DMatrix(sim.adjSens("x0"))
    adjSens_P_csim = DMatrix(sim.adjSens("p"))
    
    sol.adjSeed().set([1,0])
    sol.evaluate(0,1)
    adjSens_X0_exact = sol.adjSens(1)
    adjSens_P_exact = sol.adjSens(2)
    
    digits = 3

    self.assertAlmostEqual(adjSens_X0_int[0],adjSens_X0_exact[0],digits,"Adjoint sensitivity")
    self.assertAlmostEqual(adjSens_X0_int[1],adjSens_X0_exact[1],digits,"Adjoint sensitivity")
    self.assertAlmostEqual(adjSens_X0_sim[0],adjSens_X0_exact[0],digits,"Adjoint sensitivity")
    self.assertAlmostEqual(adjSens_X0_sim[1],adjSens_X0_exact[1],digits,"Adjoint sensitivity")
    self.assertAlmostEqual(adjSens_X0_csim[0],adjSens_X0_exact[0],digits,"Adjoint sensitivity")
    self.assertAlmostEqual(adjSens_X0_csim[1],adjSens_X0_exact[1],digits,"Adjoint sensitivity")

    self.assertAlmostEqual(adjSens_P_int[0],adjSens_P_exact[0],digits,"Adjoint sensitivity")
    self.assertAlmostEqual(adjSens_P_int[1],adjSens_P_exact[1],digits,"Adjoint sensitivity")
    self.assertAlmostEqual(adjSens_P_sim[0],adjSens_P_exact[0],digits,"Adjoint sensitivity")
    self.assertAlmostEqual(adjSens_P_sim[1],adjSens_P_exact[1],digits,"Adjoint sensitivity")
    self.assertAlmostEqual(adjSens_P_csim[0],adjSens_P_exact[0],digits,"Adjoint sensitivity")
    self.assertAlmostEqual(adjSens_P_csim[1],adjSens_P_exact[1],digits,"Adjoint sensitivity")
    
if __name__ == '__main__':
    unittest.main()

