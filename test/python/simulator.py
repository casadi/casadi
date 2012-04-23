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
    dp=ssym("dp")
    f = DAE_NUM_IN * [[]]
    f[DAE_T] = t
    f[DAE_Y] = q
    f[DAE_P] = p
    f[DAE_YDOT] = dp
    f=SXFunction(f,[q/p*t**2-dp])
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
    [qend,_]=integrator.call([q0,par,MX()])
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
    dp=ssym("dp")
    
    out = SXFunction({'NUM': DAE_NUM_IN, DAE_T: t, DAE_Y: q, DAE_P: p, DAE_YDOT: dp},[[q],[t],[p],[dp]])
    out.init()
        
    f=SXFunction({'NUM': DAE_NUM_IN, DAE_T: t, DAE_Y: q, DAE_P: p, DAE_YDOT: dp},[q/p*t**2-dp])
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
    self.checkarray(sim.output(3),DMatrix.zeros(tc.shape),"Evaluation output mismatch")
    
    f=SXFunction({'NUM': DAE_NUM_IN, DAE_T: t, DAE_Y: q, DAE_P: p},[q/p*t**2])
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
    self.checkarray(sim.output(3),DMatrix.zeros(tc.shape),"Evaluation output mismatch")
    
    out = SXFunction({'NUM': DAE_NUM_IN, DAE_T: t, DAE_Y: q, DAE_YDOT: dp},[[q],[t],[dp]])
    out.init()
    
    f=SXFunction({'NUM': DAE_NUM_IN, DAE_T: t, DAE_Y: q},[q/num['p']*t**2])
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
    self.checkarray(sim.output(2),DMatrix.zeros(tc.shape),"Evaluation output mismatch")
    
    f=SXFunction({'NUM': DAE_NUM_IN, DAE_Y: q},[-q])
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
    self.checkarray(sim.output(2),DMatrix.zeros(tc.shape),"Evaluation output mismatch")
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
    out = SXFunction({'NUM': DAE_NUM_IN, DAE_T: tv},[tv])
    
    out.init()
    sim = Simulator(self.integrator,out,t)
    sim.init()
    sim.input(0).set([num['q0']])
    sim.input(1).set([num['p']])
    sim.evaluate()

    self.checkarray(sim.output(),t,"Evaluation output mismatch")

    pv = SX("p")
    out = SXFunction({'NUM': DAE_NUM_IN, DAE_P: pv},[pv])
    
    out.init()
    sim = Simulator(self.integrator,out,t)
    sim.init()
    sim.input(0).set([num['q0']])
    sim.input(1).set([num['p']])
    sim.evaluate()

    
    self.checkarray(sim.output(),DMatrix.ones(sim.output().shape)*p,"Evaluation output mismatch")

    yv = SX("y")
    out = SXFunction({'NUM': DAE_NUM_IN, DAE_YDOT: yv},[yv])
    
    out.init()
    sim = Simulator(self.integrator,out,t)
    sim.init()
    sim.input(0).set([num['q0']])
    sim.input(1).set([num['p']])
    sim.evaluate()

    self.checkarray(sim.output(),DMatrix.zeros(sim.output().shape),"INTEGRATOR_XPF unsupported",digits=9)
    #self.checkarray(sim.output(),DMatrix((q0*t**2*exp(t**3/(3*p)))/p),"Evaluation output mismatch",digits=9)

  def test_controlsim_inputs(self):
    self.message("ControlSimulator: inputs")
    num=self.num
    tc = 0.01*DMatrix([0,8,16,24,32])
    t  = ssym("t")
    q  = ssym("q")
    p  = ssym("p")
    dp = ssym("dp")
    
    t0 = ssym("t0")
    tf_= ssym("tf")
    
    qm = ssym("qm")
    
    out = SXFunction({'NUM': CONTROL_DAE_NUM_IN, CONTROL_DAE_T: t, CONTROL_DAE_Y: q, CONTROL_DAE_P: p, CONTROL_DAE_U: [], CONTROL_DAE_U_INTERP: [],  CONTROL_DAE_YDOT: dp, CONTROL_DAE_T0: t0 , CONTROL_DAE_TF: tf_, CONTROL_DAE_Y_MAJOR: qm},[[q],[t],[p],[dp],[t0],[tf_],[qm]])
    out.init()
    
    f=SXFunction({'NUM': CONTROL_DAE_NUM_IN, CONTROL_DAE_T: t, CONTROL_DAE_Y: q, CONTROL_DAE_P: p, CONTROL_DAE_YDOT: dp, CONTROL_DAE_Y_MAJOR: qm},[q/p*t**2-dp])
    f.init()
    sim = ControlSimulator(f,out,tc)
    sim.setOption('nf',4)
    sim.setOption('integrator',CVodesIntegrator)
    sim.setOption('integrator_options', {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True})
    sim.init()
    sim.input(CONTROLSIMULATOR_X0).set([num['q0']])
    sim.input(CONTROLSIMULATOR_P).set([num['p']])
    self.assertTrue(sim.input(CONTROLSIMULATOR_U).empty())
    sim.evaluate()
    
    tf = DMatrix(sim.getMinorT())


    self.checkarray(sim.output(),num['q0']*exp(tf**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.output(1),tf,"Evaluation output mismatch")
    self.checkarray(sim.output(2),DMatrix.ones(tf.shape)*num['p'],"Evaluation output mismatch")
    self.checkarray(sim.output(3),DMatrix.zeros(tf.shape),"Evaluation output mismatch")
    self.checkarray(sim.output(4),DMatrix([0,0,0,0,0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.24]),"Evaluation output mismatch")
    self.checkarray(sim.output(5),DMatrix([0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.32,0.32,0.32,0.32,0.32]),"Evaluation output mismatch")
    self.checkarray(sim.output(6),num['q0']*exp(sim.output(4)**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    
    
    #CasadiOptions.setCatchErrorsPython(False)
    f=SXFunction({'NUM': CONTROL_DAE_NUM_IN, CONTROL_DAE_T: t, CONTROL_DAE_Y: q, CONTROL_DAE_P: p, CONTROL_DAE_Y_MAJOR: qm},[q/p*t**2])
    f.init()
    sim = ControlSimulator(f,out,tc)
    sim.setOption('nf',4)
    sim.setOption('integrator',CVodesIntegrator)
    sim.setOption('integrator_options', {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True})
    sim.init()
    sim.input(CONTROLSIMULATOR_X0).set([num['q0']])
    sim.input(CONTROLSIMULATOR_P).set([num['p']])
    sim.evaluate()
    
    tf = DMatrix(sim.getMinorT())

    self.checkarray(sim.output(),num['q0']*exp(tf**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.output(1),tf,"Evaluation output mismatch")
    self.checkarray(sim.output(2),DMatrix.ones(tf.shape)*num['p'],"Evaluation output mismatch")
    self.checkarray(sim.output(3),DMatrix.zeros(tf.shape),"Evaluation output mismatch")
    self.checkarray(sim.output(4),DMatrix([0,0,0,0,0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.24]),"Evaluation output mismatch")
    self.checkarray(sim.output(5),DMatrix([0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.32,0.32,0.32,0.32,0.32]),"Evaluation output mismatch")
    self.checkarray(sim.output(6),num['q0']*exp(sim.output(4)**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    
    out = SXFunction({'NUM': CONTROL_DAE_NUM_IN, CONTROL_DAE_T: t, CONTROL_DAE_Y: q, CONTROL_DAE_P: [], CONTROL_DAE_U: [], CONTROL_DAE_U_INTERP: [],  CONTROL_DAE_YDOT: dp, CONTROL_DAE_T0: t0 , CONTROL_DAE_TF: tf_, CONTROL_DAE_Y_MAJOR: qm},[[q],[t],[dp],[t0],[tf_],[qm]])
    out.init()
    
    f=SXFunction({'NUM': CONTROL_DAE_NUM_IN, CONTROL_DAE_T: t, CONTROL_DAE_Y: q, CONTROL_DAE_Y_MAJOR: qm},[q/num['p']*t**2])
    f.init()
    sim = ControlSimulator(f,out,tc)
    sim.setOption('nf',4)
    sim.setOption('integrator',CVodesIntegrator)
    sim.setOption('integrator_options', {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True})
    sim.init()
    sim.input(CONTROLSIMULATOR_X0).set([num['q0']])
    sim.evaluate()

    self.checkarray(sim.output(),num['q0']*exp(tf**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.output(1),tf,"Evaluation output mismatch")
    self.checkarray(sim.output(2),DMatrix.zeros(tf.shape),"Evaluation output mismatch")
    self.checkarray(sim.output(3),DMatrix([0,0,0,0,0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.24]),"Evaluation output mismatch")
    self.checkarray(sim.output(4),DMatrix([0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.32,0.32,0.32,0.32,0.32]),"Evaluation output mismatch")
    self.checkarray(sim.output(5),num['q0']*exp(sim.output(3)**3/(3*num['p'])),"Evaluation output mismatch",digits=9)
    
    f=SXFunction({'NUM': CONTROL_DAE_NUM_IN, CONTROL_DAE_Y: q, CONTROL_DAE_Y_MAJOR: qm},[-q])
    f.init()
    sim = ControlSimulator(f,out,tc)
    sim.setOption('nf',4)
    sim.setOption('integrator',CVodesIntegrator)
    sim.setOption('integrator_options', {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True})
    sim.init()
    sim.input(CONTROLSIMULATOR_X0).set([num['q0']])
    sim.evaluate()


    self.checkarray(sim.output(),num['q0']*exp(-tf),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.output(1),tf,"Evaluation output mismatch")
    self.checkarray(sim.output(2),DMatrix.zeros(tf.shape),"Evaluation output mismatch")
    self.checkarray(sim.output(3),DMatrix([0,0,0,0,0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.24]),"Evaluation output mismatch")
    self.checkarray(sim.output(4),DMatrix([0.08,0.08,0.08,0.08,0.16,0.16,0.16,0.16,0.24,0.24,0.24,0.24,0.32,0.32,0.32,0.32,0.32]),"Evaluation output mismatch")
    self.checkarray(sim.output(5),num['q0']*exp(-sim.output(3)),"Evaluation output mismatch",digits=9)

  def test_controlsim_interpolation(self):
    q  = ssym("q")
    u  = ssym("u")
    ui = ssym("ui")
    dq = ssym("dq")


    
    tc = 0.01*DMatrix([0,8,16,24,32])
    
    U = DMatrix([0,0.1,0,0.2])
    
    q0=2.3

    out = SXFunction({'NUM': CONTROL_DAE_NUM_IN, CONTROL_DAE_Y: q, CONTROL_DAE_YDOT: dq, CONTROL_DAE_U:u, CONTROL_DAE_U_INTERP: ui},[[q],[u],[ui]])
    out.init()
    
    f=SXFunction({'NUM': CONTROL_DAE_NUM_IN, CONTROL_DAE_Y: q, CONTROL_DAE_U:u},[-q])
    f.init()
    sim = ControlSimulator(f,out,tc)
    sim.setOption('nf',4)
    sim.setOption('integrator',CVodesIntegrator)
    sim.setOption('integrator_options', {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True})
    self.assertRaises(Exception,lambda : sim.init())
    
    out = SXFunction({'NUM': CONTROL_DAE_NUM_IN, CONTROL_DAE_Y: q, CONTROL_DAE_YDOT: dq, CONTROL_DAE_U:u},[[q],[u]])
    out.init()
    
    f=SXFunction({'NUM': CONTROL_DAE_NUM_IN, CONTROL_DAE_Y: q},[-q])
    f.init()
    sim = ControlSimulator(f,out,tc)
    sim.setOption('nf',4)
    sim.setOption('integrator',CVodesIntegrator)
    sim.setOption('integrator_options', {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True})
    self.assertRaises(Exception,lambda : sim.init())


    out = SXFunction({'NUM': CONTROL_DAE_NUM_IN, CONTROL_DAE_Y: q, CONTROL_DAE_YDOT: dq, CONTROL_DAE_U: u , CONTROL_DAE_U_INTERP: ui},[[q],[u],[ui]])
    out.init()
    
    f=SXFunction({'NUM': CONTROL_DAE_NUM_IN, CONTROL_DAE_Y: q,CONTROL_DAE_U:u, CONTROL_DAE_U_INTERP: ui},[-q])
    f.init()
    sim = ControlSimulator(f,out,tc)
    sim.setOption('nf',4)
    sim.setOption('integrator',CVodesIntegrator)
    sim.setOption('integrator_options', {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True})
    sim.init()
    sim.input(CONTROLSIMULATOR_X0).set([2.3])
    sim.input(CONTROLSIMULATOR_U).set(U)
    sim.evaluate()

    tf = DMatrix(sim.getMinorT())
    
    self.checkarray(sim.output(),q0*exp(-tf),"Evaluation output mismatch",digits=9)
    self.checkarray(sim.output(1),DMatrix([0,0,0,0,0.1,0.1,0.1,0.1,0,0,0,0,0.2,0.2,0.2,0.2,0.2]),"Evaluation output mismatch")
    self.checkarray(sim.output(2),DMatrix([0,0.025,0.05,0.075,0.1,0.075,0.05,0.025,0,0.05,0.1,0.15,0.2,0.2,0.2,0.2,0.2]),"Evaluation output mismatch")

    out = SXFunction({'NUM': CONTROL_DAE_NUM_IN, CONTROL_DAE_Y: q, CONTROL_DAE_YDOT: dq, CONTROL_DAE_U: u , CONTROL_DAE_U_INTERP: ui},[[q],[u],[ui]])
    out.init()
    
    f=SXFunction({'NUM': CONTROL_DAE_NUM_IN, CONTROL_DAE_Y: q,CONTROL_DAE_U:u, CONTROL_DAE_U_INTERP: ui},[-q])
    f.init()
    sim = ControlSimulator(f,out,tc)
    sim.setOption('control_endpoint',True)
    sim.setOption('nf',4)
    sim.setOption('integrator',CVodesIntegrator)
    sim.setOption('integrator_options', {"reltol":1e-15,"abstol":1e-15,"fsens_err_con": True})
    sim.init()
    sim.input(CONTROLSIMULATOR_X0).set([2.3])
    sim.input(CONTROLSIMULATOR_U).set(vertcat([U,0]))
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
    dp=ssym("dp")
    f=SXFunction({'NUM': CONTROL_DAE_NUM_IN, CONTROL_DAE_T: t, CONTROL_DAE_Y: q, CONTROL_DAE_P: p, CONTROL_DAE_YDOT: dp},[q/p*t**2-dp])
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
    t = SX("t")

    x  = SX("x") 
    dx = SX("dx") 
    v  = SX("v") 
    dv = SX("dv") 

    u = SX("u") 

    b = SX("b")
    b_ = 0.1
    k = SX("k")
    k_ = 1
    

    rhs = vertcat([v - dx, ( -  b*v - k*x) - dv ])
    f=SXFunction({'NUM': DAE_NUM_IN, DAE_T: t, DAE_YDOT: [dx,dv], DAE_Y: [x,v], DAE_P: [b,k]},[rhs])
    f.init()
    
    cf=SXFunction({'NUM': CONTROL_DAE_NUM_IN, CONTROL_DAE_T: t, CONTROL_DAE_YDOT: [dx,dv], CONTROL_DAE_Y: [x,v], CONTROL_DAE_P: [b,k]},[rhs])
    cf.init()

    x0 = SX("x0")
    dx0 = SX("dx0")


    X0 = DMatrix([1,0])
    p_ = [b_,k_]

    # algebraic solution
    sole = exp(-(b*t)/2)*((sin((sqrt(4-b**2)*t)/2)*(2*(dx0+x0*b)-x0*b))/sqrt(4-b**2)+x0*cos((sqrt(4-b**2)*t)/2))
    sol = SXFunction([[t],[x0,dx0],[b,k]],[vertcat([sole,jacobian(sole,t)])])
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

      integrator.input(INTEGRATOR_X0).set(X0)
      integrator.input(INTEGRATOR_P).set(p_)
      integrator.fwdSeed(INTEGRATOR_X0).set([1,0])
      integrator.fwdSeed(INTEGRATOR_P).set([0,0])
      integrator.evaluate(1,0)
      
      fwdSens_int = DMatrix(integrator.fwdSens(INTEGRATOR_XF))

      ts = linspace(0,50,100)

      sim=Simulator(integrator,ts)
      sim.init()
      sim.input(INTEGRATOR_X0).set(X0)
      sim.input(INTEGRATOR_P).set(p_)
      sim.fwdSeed(INTEGRATOR_X0).set([1,0])
      sim.fwdSeed(INTEGRATOR_P).set([0,0])
      sim.evaluate(1,0)

      fwdSens_sim = DMatrix(sim.fwdSens(INTEGRATOR_XF)[-1,:])
      
      csim = ControlSimulator(cf,ts)
      csim.setOption("integrator",Integrator)
      csim.setOption("integrator_options",{"fsens_abstol": 1e-12, "fsens_reltol": 1e-12, "fsens_err_con": True})
      csim.init()
      csim.input(CONTROLSIMULATOR_X0).set(X0)
      csim.input(CONTROLSIMULATOR_P).set(p_)
      csim.fwdSeed(CONTROLSIMULATOR_X0).set([1,0])
      csim.fwdSeed(CONTROLSIMULATOR_P).set([0,0])
      csim.evaluate(1,0)
      fwdSens_csim = DMatrix(sim.fwdSens(INTEGRATOR_XF)[-1,:])
      
      sol.fwdSeed(1).set([1,0])
      sol.fwdSeed(2).set([0,0])
      sol.evaluate(1,0)
      fwdSens_exact = sol.fwdSens()
      
      digits = 6
      if (Integrator is IdasIntegrator):
        digits = 2 

      self.assertAlmostEqual(fwdSens_int[0],fwdSens_exact[0],digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_int[1],fwdSens_exact[1],digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_sim[0],fwdSens_exact[0],digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_sim[1],fwdSens_exact[1],digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_csim[0],fwdSens_exact[0],digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_csim[1],fwdSens_exact[1],digits,"Forward sensitivity")
      
      integrator.fwdSeed(INTEGRATOR_X0).set([0,0])
      integrator.fwdSeed(INTEGRATOR_P).set([1,0])
      integrator.evaluate(1,0)
      
      fwdSens_int = DMatrix(integrator.fwdSens(INTEGRATOR_XF))
      
      sim.fwdSeed(INTEGRATOR_X0).set([0,0])
      sim.fwdSeed(INTEGRATOR_P).set([1,0])
      sim.evaluate(1,0)

      fwdSens_sim = DMatrix(sim.fwdSens(INTEGRATOR_XF)[-1,:])

      csim.fwdSeed(CONTROLSIMULATOR_X0).set([0,0])
      csim.fwdSeed(CONTROLSIMULATOR_P).set([1,0])
      csim.evaluate(1,0)
      fwdSens_csim = DMatrix(sim.fwdSens(INTEGRATOR_XF)[-1,:])
      
      sol.fwdSeed(1).set([0,0])
      sol.fwdSeed(2).set([1,0])
      sol.evaluate(1,0)
      fwdSens_exact = sol.fwdSens()
      
      digits = 6
      if (Integrator is IdasIntegrator):
        digits = 2 

      self.assertAlmostEqual(fwdSens_int[0],fwdSens_exact[0], digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_int[1],fwdSens_exact[1], digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_sim[0],fwdSens_exact[0], digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_sim[1],fwdSens_exact[1], digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_csim[0],fwdSens_exact[0], digits,"Forward sensitivity")
      self.assertAlmostEqual(fwdSens_csim[1],fwdSens_exact[1], digits,"Forward sensitivity")
        
    
if __name__ == '__main__':
    unittest.main()

