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



scipy_available = True
try:
	import scipy.special
	from scipy.linalg import expm
except:
	scipy_available = False

class Simulatortests(casadiTestCase):


  def setUp(self):
    # Reference solution is q0 e^((t^3-t0^3)/(3 p))
    t=SX.sym("t")
    q=SX.sym("q")
    p=SX.sym("p")
    f=SXFunction('f', daeIn(t=t,x=q,p=p),daeOut(ode=q/p*t**2))
    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    opts["fsens_err_con"] = True
    #opts["verbose"] = True
    opts["t0"] = 0
    opts["tf"] = 2.3
    integrator = Integrator("integrator", "cvodes", f, opts)
    q0   = MX.sym("q0")
    par  = MX.sym("p")
    qend = integrator({'x0':q0,'p':par})["xf"]
    qe=MXFunction('qe', [q0,par],[qend])
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
    
    t=SX.sym("t")
    q=SX.sym("q")
    p=SX.sym("p")
    
    out = SXFunction('out', daeIn(t=t, x=q, p=p),[q])
        
    f=SXFunction('f', daeIn(t=t, x=q, p=p),daeOut(ode=q/p*t**2))
    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    opts["fsens_err_con"] = True
    #opts["verbose"] = True
    opts["t0"] = 0
    opts["tf"] = 2.3
    integrator = Integrator("integrator", "cvodes", f, opts)
    sim = Simulator("sim", integrator, out, tc)
    solution = SXFunction('solution', integratorIn(x0=q, p=p),[horzcat([q*exp(t**3/(3*p)) for t in tc])])
    
    for f in [sim,solution]:
      f.setInput(0.3,"x0")
      f.setInput(0.7,"p")
    
    self.checkfunction(sim,solution,adj=False,jacobian=False,sens_der=False,evals=False,digits=6)

  def test_controlsim_full(self):
    self.message("ControlSimulator inputs")
    num = self.num
    N = 4
    tc = DMatrix(n.linspace(0,num['tend'],N))
    
    t=SX.sym("t")
    q=SX.sym("q")
    p=SX.sym("p")
    u=SX.sym("u")
    out = SXFunction('out', controldaeIn(t=t, x=q, p=p),[q])
        
    f=SXFunction('f', daeIn(t=t, x=q, p=p),daeOut(ode=q/p*t**2))
    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    opts["fsens_err_con"] = True
    #opts["verbose"] = True
    opts["t0"] = 0
    opts["tf"] = 2.3
    integrator = Integrator("integrator", "cvodes", f, opts)
    cdae = SXFunction('cdae', controldaeIn(t=t,x=q,p=p,u=u),daeOut(ode=u*q/p*t**2))    
    sim = ControlSimulator("sim", cdae,out,tc, {"integrator":"cvodes",
                                                "integrator_options":{"reltol": 1e-15,
                                                                      "abstol": 1e-15,
                                                                      "fsens_err_con": True}})
    
    U = SX.sym("U",1,N-1)
    
    result = [SX(q)]
    for i in range(N-1):
      tf = tc[i+1]
      t0 = tc[i]
      xf = lambda t,t0: exp((t**3-t0**3)/3/p*U[0,i])
      result.append(result[-1]*xf(tf,t0))
    
    solution = SXFunction('solution', controlsimulatorIn(x0=q, p=p,u=U),[horzcat(result)])
    
    for f in [sim,solution]:
      f.setInput(0.3,"x0")
      f.setInput(0.7,"p")
      f.setInput(2,"u")
      f.setInput(DMatrix(range(1,N)).T/10,"u")
    
    self.checkfunction(sim,solution,adj=False,jacobian=False,sens_der=False,evals=False,digits=6)
    
  def test_sim_inputs(self):
    self.message("Simulator inputs")
    num = self.num
    tc = DMatrix(n.linspace(0,num['tend'],100))
    
    t=SX.sym("t")
    q=SX.sym("q")
    p=SX.sym("p")
    
    out = SXFunction('out', daeIn(t=t, x=q, p=p),[q,t,p])
        
    f=SXFunction('f', daeIn(t=t, x=q, p=p),daeOut(ode=q/p*t**2))
    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    opts["fsens_err_con"] = True
    #opts["verbose"] = True
    opts["t0"] = 0
    opts["tf"] = 2.3
    integrator = Integrator("integrator", "cvodes", f, opts)
    sim = Simulator("sim", integrator, out, tc)
    sim.setInput([num['q0']],0)
    sim.setInput([num['p']],1)
    sim.evaluate()

    self.checkarray(sim.getOutput().T,num['q0']*exp(tc**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.getOutput(1).T,tc,"Evaluation output mismatch")
    self.checkarray(sim.getOutput(2).T,DMatrix.ones(tc.shape)*num['p'],"Evaluation output mismatch")
    
    f=SXFunction('f', daeIn(t=t, x=q, p=p),daeOut(ode=q/p*t**2))
    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    opts["fsens_err_con"] = True
    #opts["verbose"] = True
    opts["t0"] = 0
    opts["tf"] = 2.3
    integrator = Integrator("integrator", "cvodes", f, opts)
    sim = Simulator("sim", integrator,out,tc)
    sim.setInput([num['q0']],0)
    sim.setInput([num['p']],1)
    sim.evaluate()

    self.checkarray(sim.getOutput().T,num['q0']*exp(tc**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.getOutput(1).T,tc,"Evaluation output mismatch")
    self.checkarray(sim.getOutput(2).T,DMatrix.ones(tc.shape)*num['p'],"Evaluation output mismatch")
    
    out = SXFunction('out', daeIn(t=t, x=q),[q,t])
    
    f=SXFunction('f', daeIn(t=t, x=q),daeOut(ode=q/num['p']*t**2))
    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    opts["fsens_err_con"] = True
    #opts["verbose"] = True
    opts["t0"] = 0
    opts["tf"] = 2.3
    integrator = Integrator("integrator", "cvodes", f, opts)
    sim = Simulator("sim", integrator, out, tc)
    sim.setInput([num['q0']],0)
    sim.evaluate()

    self.checkarray(sim.getOutput().T,num['q0']*exp(tc**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.getOutput(1).T,tc,"Evaluation output mismatch")
    
    f=SXFunction('f', daeIn(x=q),daeOut(ode=-q))
    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    opts["fsens_err_con"] = True
    #opts["verbose"] = True
    opts["t0"] = 0
    opts["tf"] = 2.3
    integrator = Integrator("integrator", "cvodes", f, opts)
    sim = Simulator("sim", integrator, out, tc)
    sim.setInput([num['q0']],0)
    sim.evaluate()

    self.checkarray(sim.getOutput().T,num['q0']*exp(-tc),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.getOutput(1).T,tc,"Evaluation output mismatch",digits=9)
    self.assertTrue(sim.sparsity_out(0)==sim.sparsity_out(1))
    
  def test_sim_outputs(self):
    self.message("Simulator: outputs")
    num=self.num
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    t = n.linspace(0,num['tend'],100)
    sim = Simulator("sim", self.integrator, t)
    sim.setInput([num['q0']],0)
    sim.setInput([num['p']],1)
    sim.evaluate()

    self.checkarray(sim.getOutput().T,DMatrix(q0*exp(t**3/(3*p))),"Evaluation output mismatch",digits=9)

    tv = SX.sym("t")
    out = SXFunction('out', daeIn(t=tv),[tv])
    
    sim = Simulator('sim', self.integrator,out,t)
    sim.setInput([num['q0']],0)
    sim.setInput([num['p']],1)
    sim.evaluate()

    self.checkarray(sim.getOutput().T,t,"Evaluation output mismatch")

    pv = SX.sym("p")
    out = SXFunction('out', daeIn(p=pv),[pv])
    
    sim = Simulator('sim', self.integrator,out,t)
    sim.setInput([num['q0']],0)
    sim.setInput([num['p']],1)
    sim.evaluate()

    
    self.checkarray(sim.getOutput(),DMatrix.ones(sim.size_out(0))*p,"Evaluation output mismatch")

    #yv = SX.sym("y")
    #out = SXFunction('out', daeIn(),[yv])
    
    #sim = Simulator('sim', self.integrator,out,t)
    #sim.setInput([num['q0']],0)
    #sim.setInput([num['p']],1)
    #sim.evaluate()

    #self.checkarray(sim.getOutput(),DMatrix.zeros(sim.size_out(0)),"'xpf' unsupported",digits=9)
    #self.checkarray(sim.getOutput(),DMatrix((q0*t**2*exp(t**3/(3*p)))/p),"Evaluation output mismatch",digits=9)

  def test_controlsim_inputs(self):
    self.message("ControlSimulator: inputs")
    num=self.num
    tc = 0.01*DMatrix([0,8,16,24,32])
    t  = SX.sym("t")
    q  = SX.sym("q")
    p  = SX.sym("p")
    
    t0 = SX.sym("t0")
    tf_= SX.sym("tf")
    
    qm = SX.sym("qm")
    
    out = SXFunction('out', controldaeIn(t=t, x=q, p=p, t0=t0, tf=tf_, x_major=qm),[q,t,p,t0,tf_,qm])
    
    f=SXFunction('f', controldaeIn(t=t, x=q, p=p, x_major=qm),[q/p*t**2])
    opts = {}
    opts['nf'] = 4
    opts['integrator'] = "cvodes"
    opts['integrator_options'] = {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True}
    sim = ControlSimulator("sim", f, out, tc, opts)
    sim.setInput([num['q0']],"x0")
    sim.setInput([num['p']],"p")
    self.assertTrue(sim.sparsity_in("u").isempty())
    sim.evaluate()
    
    tf = DMatrix(sim.getMinorT())


    self.checkarray(sim.getOutput(),num['q0']*exp(tf**3/(3*num['p'])).T,"Evaluation output mismatch",digits=9)
    self.checkarray(sim.getOutput(1),tf.T,"Evaluation output mismatch")
    self.checkarray(sim.getOutput(2),DMatrix.ones(tf.shape).T*num['p'],"Evaluation output mismatch")
    self.checkarray(sim.getOutput(3),DMatrix([0,0,0,0,0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.24]).T,"Evaluation output mismatch")
    self.checkarray(sim.getOutput(4),DMatrix([0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.32,0.32,0.32,0.32,0.32]).T,"Evaluation output mismatch")
    self.checkarray(sim.getOutput(5),num['q0']*exp(sim.getOutput(3)**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    
    
    #CasadiOptions.setCatchErrorsPython(False)
    f=SXFunction('f', controldaeIn(t=t, x=q, p=p, x_major=qm),[q/p*t**2])
    opts = {}
    opts['nf'] = 4
    opts['integrator'] = "cvodes"
    opts['integrator_options'] = {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True}
    sim = ControlSimulator("sim", f, out, tc, opts)
    sim.setInput([num['q0']],"x0")
    sim.setInput([num['p']],"p")
    sim.evaluate()
    
    tf = DMatrix(sim.getMinorT())

    self.checkarray(sim.getOutput().T,num['q0']*exp(tf**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.getOutput(1).T,tf,"Evaluation output mismatch")
    self.checkarray(sim.getOutput(2).T,DMatrix.ones(tf.shape)*num['p'],"Evaluation output mismatch")
    self.checkarray(sim.getOutput(3).T,DMatrix([0,0,0,0,0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.24]),"Evaluation output mismatch")
    self.checkarray(sim.getOutput(4).T,DMatrix([0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.32,0.32,0.32,0.32,0.32]),"Evaluation output mismatch")
    self.checkarray(sim.getOutput(5),num['q0']*exp(sim.getOutput(3)**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    
    out = SXFunction('out', controldaeIn(t=t, x=q, t0=t0, tf=tf_, x_major=qm),[q,t,t0,tf_,qm])
    
    f=SXFunction('f', controldaeIn(t=t, x=q, x_major=qm),[q/num['p']*t**2])
    opts = {}
    opts['nf'] = 4
    opts['integrator'] = "cvodes"
    opts['integrator_options'] = {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True}
    sim = ControlSimulator("sim", f,out, tc, opts)
    sim.setInput([num['q0']],"x0")
    sim.evaluate()

    self.checkarray(sim.getOutput().T,num['q0']*exp(tf**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.getOutput(1).T,tf,"Evaluation output mismatch")
    self.checkarray(sim.getOutput(2).T,DMatrix([0,0,0,0,0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.24]),"Evaluation output mismatch")
    self.checkarray(sim.getOutput(3).T,DMatrix([0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.32,0.32,0.32,0.32,0.32]),"Evaluation output mismatch")
    self.checkarray(sim.getOutput(4),num['q0']*exp(sim.getOutput(2)**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    
    f=SXFunction('f', controldaeIn(x=q, x_major=qm),[-q])
    opts = {}
    opts['nf'] = 4
    opts['integrator']  = "cvodes"
    opts['integrator_options'] = {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True}
    sim = ControlSimulator("sim", f, out, tc, opts)
    sim.setInput([num['q0']],"x0")
    sim.evaluate()


    self.checkarray(sim.getOutput().T,num['q0']*exp(-tf),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.getOutput(1).T,tf,"Evaluation output mismatch")
    self.checkarray(sim.getOutput(2).T,DMatrix([0,0,0,0,0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.24]),"Evaluation output mismatch")
    self.checkarray(sim.getOutput(3).T,DMatrix([0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.32,0.32,0.32,0.32,0.32]),"Evaluation output mismatch")
    self.checkarray(sim.getOutput(4),num['q0']*exp(-sim.getOutput(2)),"Evaluation output mismatch",digits=9)

  def test_controlsim_interpolation(self):
    q  = SX.sym("q")
    u  = SX.sym("u")
    ui = SX.sym("ui")


    
    tc = 0.01*DMatrix([0,8,16,24,32])
    
    U = DMatrix([0,0.1,0,0.2]).T
    
    q0=2.3

    out = SXFunction('out', controldaeIn(x=q, u=u, u_interp=ui),[q,u,ui])
    
    f=SXFunction('f', controldaeIn(x=q, u=u),[-q])
    opts = {}
    opts['nf'] = 4
    opts['integrator'] = "cvodes"
    opts['integrator_options'] = {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True}
    self.assertRaises(Exception,lambda : ControlSimulator("sim", f, out, tc, opts))
    
    out = SXFunction('out', controldaeIn(x=q, u=u),[q,u])
    
    f=SXFunction('f', controldaeIn(x=q),[-q])
    opts = {}
    opts['nf'] = 4
    opts['integrator'] = "cvodes"
    opts['integrator_options'] = {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True}
    self.assertRaises(Exception,lambda : ControlSimulator("sim", f, out, tc, opts))

    out = SXFunction('out', controldaeIn(x=q, u=u, u_interp=ui),[q,u,ui])
    
    f=SXFunction('f', controldaeIn(x=q, u=u, u_interp=ui),[-q])
    opts = {}
    opts['nf'] = 4
    opts['integrator'] = "cvodes"
    opts['integrator_options'] = {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True}
    sim = ControlSimulator("sim", f, out, tc, opts)
    sim.setInput([2.3],"x0")
    sim.setInput(U,"u")
    sim.evaluate()

    tf = DMatrix(sim.getMinorT())
    
    self.checkarray(sim.getOutput().T,q0*exp(-tf),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.getOutput(1).T,DMatrix([0,0,0,0,0.1,0.1,0.1,0.1,0,0,0,0,0.2,0.2,0.2,0.2,0.2]),"Evaluation output mismatch")
    self.checkarray(sim.getOutput(2).T,DMatrix([0,0.025,0.05,0.075,0.1,0.075,0.05,0.025,0,0.05,0.1,0.15,0.2,0.2,0.2,0.2,0.2]),"Evaluation output mismatch")

    out = SXFunction('out', controldaeIn(x=q, u=u, u_interp=ui),[q,u,ui])
    
    f=SXFunction('f', controldaeIn(x=q,u=u, u_interp=ui),[-q])
    opts = {}
    opts['control_endpoint'] = True
    opts['nf'] = 4
    opts['integrator'] = "cvodes"
    opts['integrator_options'] = {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True}
    sim = ControlSimulator("sim", f, out, tc, opts)
    sim.setInput([2.3],"x0")
    sim.setInput(horzcat([U,0]),"u")
    sim.evaluate()

    tf = DMatrix(sim.getMinorT())
    
    self.checkarray(sim.getOutput().T,q0*exp(-tf),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.getOutput(1).T,DMatrix([0,0,0,0,0.1,0.1,0.1,0.1,0,0,0,0,0.2,0.2,0.2,0.2,0.2]),"Evaluation output mismatch")
    self.checkarray(sim.getOutput(2).T,DMatrix([0,0.025,0.05,0.075,0.1,0.075,0.05,0.025,0,0.05,0.1,0.15,0.2,0.15,0.1,0.05,0]),"Evaluation output mismatch")
    #self.checkarray(sim.getMinorU(),DMatrix([]),"Evaluation output mismatch")
    

  def test_controlsim_outputs(self):
    self.message("CVodes integration: outputs")
    num=self.num
    tc = 0.01*DMatrix([0,8,16,24,32])
    t=SX.sym("t")
    q=SX.sym("q")
    p=SX.sym("p")
    f=SXFunction('f', controldaeIn(t=t, x=q, p=p),[q/p*t**2])
    opts = {}
    opts['nf'] = 4
    opts['integrator'] = "cvodes"
    opts['integrator_options'] = {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True}
    sim = ControlSimulator("sim", f, tc, opts)
    sim.setInput([num['q0']],0)
    sim.setInput([num['p']],1)
    sim.evaluate()
    
    tf = DMatrix(sim.getMinorT())

    self.checkarray(sim.getOutput().T,num['q0']*exp(tf**3/(3*num['p'])),"Evaluation output mismatch",digits=9)


  def test_simulator_time_offset(self):
    self.message("CVodes integration: simulator time offset")
    num=self.num
    t = n.linspace(0.7,num['tend'],100)
    sim = Simulator("sim", self.integrator, t)
    sim.setInput([num['q0']],0)
    sim.setInput([num['p']],1)
    sim.evaluate()

    tend=num['tend']
    q0=num['q0']
    p=num['p']

    self.assertAlmostEqual(sim.getOutput()[0,-1],q0*exp((tend**3-0.7**3)/(3*p)),9,"Evaluation output mismatch")
            
if __name__ == '__main__':
    unittest.main()

