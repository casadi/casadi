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
    f=SX.fun('f', daeIn(t=t,x=q,p=p),daeOut(ode=q/p*t**2))
    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    opts["fsens_err_con"] = True
    #opts["verbose"] = True
    opts["t0"] = 0
    opts["tf"] = 2.3
    integrator = Function.integrator("integrator", "cvodes", f, opts)
    q0   = MX.sym("q0")
    par  = MX.sym("p")
    qend = integrator({'x0':q0,'p':par})["xf"]
    qe=MX.fun('qe', [q0,par],[qend])
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
    
    out = SX.fun('out', daeIn(t=t, x=q, p=p),[q])
        
    f=SX.fun('f', daeIn(t=t, x=q, p=p),daeOut(ode=q/p*t**2))
    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    opts["fsens_err_con"] = True
    #opts["verbose"] = True
    opts["t0"] = 0
    opts["tf"] = 2.3
    integrator = Function.integrator("integrator", "cvodes", f, opts)
    sim = Simulator("sim", integrator, out, tc)
    solution = SX.fun('solution', integratorIn(x0=q, p=p),[horzcat([q*exp(t**3/(3*p)) for t in tc])])
    
    for f in [sim,solution]:
      f.setInput(0.3,"x0")
      f.setInput(0.7,"p")
    
    self.checkfunction(sim,solution,adj=False,jacobian=False,sens_der=False,evals=False,digits=6)

  def test_sim_inputs(self):
    self.message("Simulator inputs")
    num = self.num
    tc = DMatrix(n.linspace(0,num['tend'],100))
    
    t=SX.sym("t")
    q=SX.sym("q")
    p=SX.sym("p")
    
    out = SX.fun('out', daeIn(t=t, x=q, p=p),[q,t,p])
        
    f=SX.fun('f', daeIn(t=t, x=q, p=p),daeOut(ode=q/p*t**2))
    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    opts["fsens_err_con"] = True
    #opts["verbose"] = True
    opts["t0"] = 0
    opts["tf"] = 2.3
    integrator = Function.integrator("integrator", "cvodes", f, opts)
    sim = Simulator("sim", integrator, out, tc)
    sim.setInput([num['q0']],0)
    sim.setInput([num['p']],1)
    sim.evaluate()

    self.checkarray(sim.getOutput().T,num['q0']*exp(tc**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.getOutput(1).T,tc,"Evaluation output mismatch")
    self.checkarray(sim.getOutput(2).T,DMatrix.ones(tc.shape)*num['p'],"Evaluation output mismatch")
    
    f=SX.fun('f', daeIn(t=t, x=q, p=p),daeOut(ode=q/p*t**2))
    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    opts["fsens_err_con"] = True
    #opts["verbose"] = True
    opts["t0"] = 0
    opts["tf"] = 2.3
    integrator = Function.integrator("integrator", "cvodes", f, opts)
    sim = Simulator("sim", integrator,out,tc)
    sim.setInput([num['q0']],0)
    sim.setInput([num['p']],1)
    sim.evaluate()

    self.checkarray(sim.getOutput().T,num['q0']*exp(tc**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.getOutput(1).T,tc,"Evaluation output mismatch")
    self.checkarray(sim.getOutput(2).T,DMatrix.ones(tc.shape)*num['p'],"Evaluation output mismatch")
    
    out = SX.fun('out', daeIn(t=t, x=q),[q,t])
    
    f=SX.fun('f', daeIn(t=t, x=q),daeOut(ode=q/num['p']*t**2))
    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    opts["fsens_err_con"] = True
    #opts["verbose"] = True
    opts["t0"] = 0
    opts["tf"] = 2.3
    integrator = Function.integrator("integrator", "cvodes", f, opts)
    sim = Simulator("sim", integrator, out, tc)
    sim.setInput([num['q0']],0)
    sim.evaluate()

    self.checkarray(sim.getOutput().T,num['q0']*exp(tc**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.getOutput(1).T,tc,"Evaluation output mismatch")
    
    f=SX.fun('f', daeIn(x=q),daeOut(ode=-q))
    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    opts["fsens_err_con"] = True
    #opts["verbose"] = True
    opts["t0"] = 0
    opts["tf"] = 2.3
    integrator = Function.integrator("integrator", "cvodes", f, opts)
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
    out = SX.fun('out', daeIn(t=tv),[tv])
    
    sim = Simulator('sim', self.integrator,out,t)
    sim.setInput([num['q0']],0)
    sim.setInput([num['p']],1)
    sim.evaluate()

    self.checkarray(sim.getOutput().T,t,"Evaluation output mismatch")

    pv = SX.sym("p")
    out = SX.fun('out', daeIn(p=pv),[pv])
    
    sim = Simulator('sim', self.integrator,out,t)
    sim.setInput([num['q0']],0)
    sim.setInput([num['p']],1)
    sim.evaluate()

    
    self.checkarray(sim.getOutput(),DMatrix.ones(sim.size_out(0))*p,"Evaluation output mismatch")

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

