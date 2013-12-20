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

  def test_sim_full(self):
    self.message("Simulator inputs")
    num = self.num
    N = 4
    tc = DMatrix(n.linspace(0,num['tend'],N))
    
    t=ssym("t")
    q=ssym("q")
    p=ssym("p")
    
    out = SXFunction(daeIn(t=t, x=q, p=p),[q])
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
    
    solution = SXFunction(integratorIn(x0=q, p=p),[vertcat([q*exp(t**3/(3*p)) for t in tc])])
    solution.init()
    
    for f in [sim,solution]:
      f.setInput(0.3,"x0")
      f.setInput(0.7,"p")
    
    self.checkfx(sim,solution,adj=False,sens_der=False,evals=False,digits=6)

  def test_controlsim_full(self):
    self.message("ControlSimulator inputs")
    num = self.num
    N = 4
    tc = DMatrix(n.linspace(0,num['tend'],N))
    
    t=ssym("t")
    q=ssym("q")
    p=ssym("p")
    u=ssym("u")
    out = SXFunction(controldaeIn(t=t, x=q, p=p),[q])
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
    
    cdae = SXFunction(controldaeIn(t=t,x=q,p=p,u=u),daeOut(ode=u*q/p*t**2))
    cdae.init()
    
    sim = ControlSimulator(cdae,out,tc)
    sim.setOption("integrator",CVodesIntegrator)
    sim.setOption("integrator_options",{"reltol": 1e-15, "abstol": 1e-15, "fsens_err_con": True})
    sim.init()
    
    U = ssym("U",N-1)
    
    result = SXMatrix(q)
    for i in range(N-1):
      tf = tc[i+1]
      t0 = tc[i]
      xf = lambda t,t0: exp((t**3-t0**3)/3/p*U[i])
      result.append(result[-1]*xf(tf,t0))
    
    solution = SXFunction(controlsimulatorIn(x0=q, p=p,u=U),[result])
    solution.init()
    
    for f in [sim,solution]:
      f.setInput(0.3,"x0")
      f.setInput(0.7,"p")
      f.setInput(2,"u")
      f.setInput(DMatrix(range(1,N))/10,"u")
    
    self.checkfx(sim,solution,adj=False,sens_der=False,evals=False,digits=6)
    
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
    sim.setInput([num['q0']],0)
    sim.setInput([num['p']],1)
    sim.evaluate()

    self.checkarray(sim.getOutput(),num['q0']*exp(tc**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.getOutput(1),tc,"Evaluation output mismatch")
    self.checkarray(sim.getOutput(2),DMatrix.ones(tc.shape)*num['p'],"Evaluation output mismatch")
    
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
    sim.setInput([num['q0']],0)
    sim.setInput([num['p']],1)
    sim.evaluate()

    self.checkarray(sim.getOutput(),num['q0']*exp(tc**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.getOutput(1),tc,"Evaluation output mismatch")
    self.checkarray(sim.getOutput(2),DMatrix.ones(tc.shape)*num['p'],"Evaluation output mismatch")
    
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
    sim.setInput([num['q0']],0)
    sim.evaluate()

    self.checkarray(sim.getOutput(),num['q0']*exp(tc**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.getOutput(1),tc,"Evaluation output mismatch")
    
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
    sim.setInput([num['q0']],0)
    sim.evaluate()

    self.checkarray(sim.getOutput(),num['q0']*exp(-tc),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.getOutput(1),tc,"Evaluation output mismatch",digits=9)
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
    sim.setInput([num['q0']],0)
    sim.setInput([num['p']],1)
    sim.evaluate()

    self.checkarray(sim.getOutput(),DMatrix(q0*exp(t**3/(3*p))),"Evaluation output mismatch",digits=9)

    tv = SX("t")
    out = SXFunction(daeIn(t=tv),[tv])
    
    out.init()
    sim = Simulator(self.integrator,out,t)
    sim.init()
    sim.setInput([num['q0']],0)
    sim.setInput([num['p']],1)
    sim.evaluate()

    self.checkarray(sim.getOutput(),t,"Evaluation output mismatch")

    pv = SX("p")
    out = SXFunction(daeIn(p=pv),[pv])
    
    out.init()
    sim = Simulator(self.integrator,out,t)
    sim.init()
    sim.setInput([num['q0']],0)
    sim.setInput([num['p']],1)
    sim.evaluate()

    
    self.checkarray(sim.getOutput(),DMatrix.ones(sim.output().shape)*p,"Evaluation output mismatch")

    #yv = SX("y")
    #out = SXFunction(daeIn(),[yv])
    
    #out.init()
    #sim = Simulator(self.integrator,out,t)
    #sim.init()
    #sim.setInput([num['q0']],0)
    #sim.setInput([num['p']],1)
    #sim.evaluate()

    #self.checkarray(sim.getOutput(),DMatrix.zeros(sim.output().shape),"INTEGRATOR_XPF unsupported",digits=9)
    #self.checkarray(sim.getOutput(),DMatrix((q0*t**2*exp(t**3/(3*p)))/p),"Evaluation output mismatch",digits=9)

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
    sim.setInput([num['q0']],"x0")
    sim.setInput([num['p']],"p")
    self.assertTrue(sim.input("u").empty())
    sim.evaluate()
    
    tf = DMatrix(sim.getMinorT())


    self.checkarray(sim.getOutput(),num['q0']*exp(tf**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.getOutput(1),tf,"Evaluation output mismatch")
    self.checkarray(sim.getOutput(2),DMatrix.ones(tf.shape)*num['p'],"Evaluation output mismatch")
    self.checkarray(sim.getOutput(3),DMatrix([0,0,0,0,0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.24]),"Evaluation output mismatch")
    self.checkarray(sim.getOutput(4),DMatrix([0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.32,0.32,0.32,0.32,0.32]),"Evaluation output mismatch")
    self.checkarray(sim.getOutput(5),num['q0']*exp(sim.getOutput(3)**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    
    
    #CasadiOptions.setCatchErrorsPython(False)
    f=SXFunction(controldaeIn(t=t, x=q, p=p, x_major=qm),[q/p*t**2])
    f.init()
    sim = ControlSimulator(f,out,tc)
    sim.setOption('nf',4)
    sim.setOption('integrator',CVodesIntegrator)
    sim.setOption('integrator_options', {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True})
    sim.init()
    sim.setInput([num['q0']],"x0")
    sim.setInput([num['p']],"p")
    sim.evaluate()
    
    tf = DMatrix(sim.getMinorT())

    self.checkarray(sim.getOutput(),num['q0']*exp(tf**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.getOutput(1),tf,"Evaluation output mismatch")
    self.checkarray(sim.getOutput(2),DMatrix.ones(tf.shape)*num['p'],"Evaluation output mismatch")
    self.checkarray(sim.getOutput(3),DMatrix([0,0,0,0,0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.24]),"Evaluation output mismatch")
    self.checkarray(sim.getOutput(4),DMatrix([0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.32,0.32,0.32,0.32,0.32]),"Evaluation output mismatch")
    self.checkarray(sim.getOutput(5),num['q0']*exp(sim.getOutput(3)**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    
    out = SXFunction(controldaeIn(t=t, x=q, t0=t0, tf=tf_, x_major=qm),[q,t,t0,tf_,qm])
    out.init()
    
    f=SXFunction(controldaeIn(t=t, x=q, x_major=qm),[q/num['p']*t**2])
    f.init()
    sim = ControlSimulator(f,out,tc)
    sim.setOption('nf',4)
    sim.setOption('integrator',CVodesIntegrator)
    sim.setOption('integrator_options', {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True})
    sim.init()
    sim.setInput([num['q0']],"x0")
    sim.evaluate()

    self.checkarray(sim.getOutput(),num['q0']*exp(tf**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.getOutput(1),tf,"Evaluation output mismatch")
    self.checkarray(sim.getOutput(2),DMatrix([0,0,0,0,0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.24]),"Evaluation output mismatch")
    self.checkarray(sim.getOutput(3),DMatrix([0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.32,0.32,0.32,0.32,0.32]),"Evaluation output mismatch")
    self.checkarray(sim.getOutput(4),num['q0']*exp(sim.getOutput(2)**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    
    f=SXFunction(controldaeIn(x=q, x_major=qm),[-q])
    f.init()
    sim = ControlSimulator(f,out,tc)
    sim.setOption('nf',4)
    sim.setOption('integrator',CVodesIntegrator)
    sim.setOption('integrator_options', {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True})
    sim.init()
    sim.setInput([num['q0']],"x0")
    sim.evaluate()


    self.checkarray(sim.getOutput(),num['q0']*exp(-tf),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.getOutput(1),tf,"Evaluation output mismatch")
    self.checkarray(sim.getOutput(2),DMatrix([0,0,0,0,0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.24]),"Evaluation output mismatch")
    self.checkarray(sim.getOutput(3),DMatrix([0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.32,0.32,0.32,0.32,0.32]),"Evaluation output mismatch")
    self.checkarray(sim.getOutput(4),num['q0']*exp(-sim.getOutput(2)),"Evaluation output mismatch",digits=9)

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
    sim.setInput([2.3],"x0")
    sim.setInput(U,"u")
    sim.evaluate()

    tf = DMatrix(sim.getMinorT())
    
    self.checkarray(sim.getOutput(),q0*exp(-tf),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.getOutput(1),DMatrix([0,0,0,0,0.1,0.1,0.1,0.1,0,0,0,0,0.2,0.2,0.2,0.2,0.2]),"Evaluation output mismatch")
    self.checkarray(sim.getOutput(2),DMatrix([0,0.025,0.05,0.075,0.1,0.075,0.05,0.025,0,0.05,0.1,0.15,0.2,0.2,0.2,0.2,0.2]),"Evaluation output mismatch")

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
    sim.setInput([2.3],"x0")
    sim.setInput(vertcat([U,0]),"u")
    sim.evaluate()

    tf = DMatrix(sim.getMinorT())
    
    self.checkarray(sim.getOutput(),q0*exp(-tf),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.getOutput(1),DMatrix([0,0,0,0,0.1,0.1,0.1,0.1,0,0,0,0,0.2,0.2,0.2,0.2,0.2]),"Evaluation output mismatch")
    self.checkarray(sim.getOutput(2),DMatrix([0,0.025,0.05,0.075,0.1,0.075,0.05,0.025,0,0.05,0.1,0.15,0.2,0.15,0.1,0.05,0]),"Evaluation output mismatch")

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
    sim.setInput([num['q0']],0)
    sim.setInput([num['p']],1)
    sim.evaluate()
    
    tf = DMatrix(sim.getMinorT())

    self.checkarray(sim.getOutput(),num['q0']*exp(tf**3/(3*num['p'])),"Evaluation output mismatch",digits=9)


  def test_simulator_time_offset(self):
    self.message("CVodes integration: simulator time offset")
    num=self.num
    t = n.linspace(0.7,num['tend'],100)
    sim = Simulator(self.integrator,t)
    sim.init()
    sim.setInput([num['q0']],0)
    sim.setInput([num['p']],1)
    sim.evaluate()

    tend=num['tend']
    q0=num['q0']
    p=num['p']

    self.assertAlmostEqual(sim.getOutput()[-1],q0*exp((tend**3-0.7**3)/(3*p)),9,"Evaluation output mismatch")
    
    
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
    sol.setInput(50,0)
    sol.setInput(X0,1)
    sol.setInput(p_,2)


    for Integrator in [CVodesIntegrator, IdasIntegrator]:
      integrator = Integrator(f)
      #integrator.setOption("verbose",True)
      #integrator.setOption("monitor",["integrate"])
      integrator.setOption("fsens_abstol",1e-12)
      integrator.setOption("fsens_reltol",1e-12)
      integrator.setOption("fsens_err_con", True)
      integrator.setOption("tf",50)
      integrator.init()

      integrator.setInput(X0,"x0")
      integrator.setInput(p_,"p")
      integrator.setFwdSeed([1,0],"x0")
      integrator.setFwdSeed([0,0],"p")
      integrator.evaluate(1,0)
      
      fwdSens_int = DMatrix(integrator.getFwdSens("xf"))

      ts = linspace(0,50,100)

      sim=Simulator(integrator,ts)
      sim.init()
      sim.setInput(X0,"x0")
      sim.setInput(p_,"p")
      sim.setFwdSeed([1,0],"x0")
      sim.setFwdSeed([0,0],"p")
      sim.evaluate(1,0)

      fwdSens_sim = DMatrix(sim.getFwdSens("xf")[-1,:])
      
      csim = ControlSimulator(cf,ts)
      csim.setOption("integrator",Integrator)
      csim.setOption("integrator_options",{"fsens_abstol": 1e-12, "fsens_reltol": 1e-12, "fsens_err_con": True})
      csim.init()
      csim.setInput(X0,"x0")
      csim.setInput(p_,"p")
      csim.setFwdSeed([1,0],"x0")
      csim.setFwdSeed([0,0],"p")
      csim.evaluate(1,0)
      fwdSens_csim = DMatrix(csim.getFwdSens("xf")[-1,:])
      
      sol.setFwdSeed([1,0],1)
      sol.setFwdSeed([0,0],2)
      sol.evaluate(1,0)
      fwdSens_exact = sol.getFwdSens()
      
      #digits = 6
      #if (Integrator is IdasIntegrator):
      digits = 2 # Joel: No reason to treat Idas differently from CVodes

      self.assertAlmostEqual(fwdSens_int[0],fwdSens_exact[0],digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_int[1],fwdSens_exact[1],digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_sim[0],fwdSens_exact[0],digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_sim[1],fwdSens_exact[1],digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_csim[0],fwdSens_exact[0],digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_csim[1],fwdSens_exact[1],digits,"Forward sensitivity")
      
      integrator.setFwdSeed([0,0],"x0")
      integrator.setFwdSeed([1,0],"p")
      integrator.evaluate(1,0)
      
      fwdSens_int = DMatrix(integrator.getFwdSens("xf"))
      
      sim.setFwdSeed([0,0],"x0")
      sim.setFwdSeed([1,0],"p")
      sim.evaluate(1,0)

      fwdSens_sim = DMatrix(sim.getFwdSens("xf")[-1,:])

      csim.setFwdSeed([0,0],"x0")
      csim.setFwdSeed([1,0],"p")
      csim.evaluate(1,0)
      fwdSens_csim = DMatrix(sim.getFwdSens("xf")[-1,:])
      
      sol.setFwdSeed([0,0],1)
      sol.setFwdSeed([1,0],2)
      sol.evaluate(1,0)
      fwdSens_exact = sol.getFwdSens()
      
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
    sol.setInput(Te,0)
    sol.setInput(X0,1)
    sol.setInput(p_,2)

    Integrator = CVodesIntegrator
  
    integrator = Integrator(f)
    #integrator.setOption("verbose",True)
    #integrator.setOption("monitor",["integrate"])
    integrator.setOption("abstolB",1e-12)
    integrator.setOption("reltolB",1e-12)
    integrator.setOption("fsens_err_con", True)
    integrator.setOption("tf",Te)
    integrator.init()

    integrator.setInput(X0,"x0")
    integrator.setInput(p_,"p")
    integrator.setAdjSeed([1,0],"xf")
    integrator.evaluate(0,1)
    
    adjSens_X0_int = DMatrix(integrator.getAdjSens("x0"))
    adjSens_P_int = DMatrix(integrator.getAdjSens("p"))

    N = 100
    ts = linspace(0,Te,N)

    sim=Simulator(integrator,ts)
    sim.init()
    sim.setInput(X0,"x0")
    sim.setInput(p_,"p")
    sim.setAdjSeed([0,0]*(N-1) + [1,0],0)
    sim.evaluate(0,1)

    adjSens_X0_sim = DMatrix(sim.getAdjSens("x0"))
    adjSens_P_sim = DMatrix(sim.getAdjSens("p"))

    
    csim = ControlSimulator(cf,ts)
    csim.setOption("integrator",Integrator)
    csim.setOption("integrator_options",{"abstolB": 1e-12, "reltolB": 1e-12, "fsens_err_con": True, "steps_per_checkpoint":1000})
    csim.init()
    csim.setInput(X0,"x0")
    csim.setInput(p_,"p")
    csim.setAdjSeed([0,0]*(N-1) + [1,0],0)
    csim.evaluate(0,1)
    adjSens_X0_csim = DMatrix(sim.getAdjSens("x0"))
    adjSens_P_csim = DMatrix(sim.getAdjSens("p"))
    
    sol.setAdjSeed([1,0])
    sol.evaluate(0,1)
    adjSens_X0_exact = sol.getAdjSens(1)
    adjSens_P_exact = sol.getAdjSens(2)
    
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

