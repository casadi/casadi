from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *

class Integrationtests(casadiTestCase):

  def setUp(self):
    t=symbolic("t")
    q=symbolic("q")
    p=symbolic("p")
    f = ODE_NUM_IN * [[]]
    f[ODE_T] = t
    f[ODE_Y] = q
    f[ODE_P] = p
    f=SXFunction(f,[q/p*t**2])
    f.init()
    integrator = CVodesIntegrator(f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    integrator.setOption("verbose",True)
    integrator.init()
    t0   = MX(0)
    tend = MX("tend")
    q0   = MX("q0")
    par  = MX("p")
    qend=integrator([t0,tend,q0,par,MX(),MX()])
    qe=MXFunction([tend,q0,par],[qend])
    qe.init()
    self.integrator = integrator
    self.qe=qe
    self.num={'tend':2.3,'q0':7.1,'p':2}
    pass
    
    
  
  def test_eval(self):
    self.message('IPOPT integration: evaluation')
    num=self.num
    qe=self.qe
    qe.input(0).set([num['tend']])
    qe.input(1).set([num['q0']])
    qe.input(2).set([num['p']])
    qe.evaluate()

    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(qe.output()[0],q0*exp(tend**3/(3*p)),9,"Evaluation output mismatch")

  def test_jac0(self):
    return # jacobian wrt end time is not supported
    self.message('IPOPT integration: jacobian to end time')
    num=self.num
    J=self.qe.jacobian(0)
    J.init()
    J.input(0).set([num['tend']])
    J.input(1).set([num['q0']])
    J.input(2).set([num['p']])
    J.evaluate()
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    print J.output()[0]
    self.assertAlmostEqual(J.output()[0],(q0*tend**2*exp(tend**3/(3*p)))/p,9,"Evaluation output mismatch")
    
  def test_jac1(self):
    self.message('IPOPT integration: jacobian to q0')
    num=self.num
    J=self.qe.jacobian(1)
    J.init()
    J.input(0).set([num['tend']])
    J.input(1).set([num['q0']])
    J.input(2).set([num['p']])
    J.evaluate()
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(J.output()[0],exp(tend**3/(3*p)),9,"Evaluation output mismatch")
    
  def test_jac2(self):
    self.message('IPOPT integration: jacobian to p')
    num=self.num
    J=self.qe.jacobian(2)
    J.init()
    J.input(0).set([num['tend']])
    J.input(1).set([num['q0']])
    J.input(2).set([num['p']])
    J.evaluate()
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(J.output()[0],-(q0*tend**3*exp(tend**3/(3*p)))/(3*p**2),9,"Evaluation output mismatch")
    
  def test_hess(self):
    self.message('IPOPT integration: hessian to p: fwd-over-adjoint on integrator')
    num=self.num
    J=self.integrator.jacobian(INTEGRATOR_P,INTEGRATOR_XF)
    J.setOption("number_of_fwd_dir",0)
    J.setOption("number_of_adj_dir",1)
    J.init()
    J.input(INTEGRATOR_TF).set([num['tend']])
    J.input(INTEGRATOR_X0).set([num['q0']])
    J.input(INTEGRATOR_P).set([num['p']])
    J.adjSeed(0).set([1])
    # Evaluate
    J.evaluate(0,1)
      
    tend=num['tend']
    q0=num['q0']
    p=num['p']

    self.assertAlmostEqual(J.adjSens(INTEGRATOR_P)[0],(q0*tend**6*exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*exp(tend**3/(3*p)))/(3*p**3),9,"Evaluation output mismatch")
    
  def test_hess2(self):
    return # known issue
    t=symbolic("t")
    q=symbolic("q")
    p=symbolic("p")
    f = ODE_NUM_IN * [[]]
    f[ODE_T] = t
    f[ODE_Y] = q
    f[ODE_P] = p
    f=SXFunction(f,[q/p*t**2])
    f.init()
    integrator = CVodesIntegrator(f)
    integrator.init()
    t0   = MX(0)
    tend = MX("tend")
    q0   = MX("q0")
    par  = MX("p")
    qend=integrator([t0,tend,q0,par,MX(),MX()])
    qe=MXFunction([tend,q0,par],[qend])
    qe.init()
    J=qe.jacobian(2)
    J.init()
    J.input(0).set([2.3])
    J.input(1).set([7.1])
    J.input(2).set([2])
    J.adjSeed(0).set([1])
    J.evaluate(0,1)
    num=self.num
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(J.adjSens()[0],(q0*tend**6*exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*exp(tend**3/(3*p)))/(3*p**3),9,"Evaluation output mismatch")

  def test_hess3(self):
    self.message('IPOPT integration: hessian to p: Jacobian of integrator.jacobian')
    num=self.num
    J=self.integrator.jacobian(INTEGRATOR_P,INTEGRATOR_XF)
    J.init()
    H=Jacobian(J,INTEGRATOR_P)
    H.setOption("ad_mode","adjoint")
    H.init()
    H.input(INTEGRATOR_TF).set([num['tend']])
    H.input(INTEGRATOR_X0).set([num['q0']])
    H.input(INTEGRATOR_P).set([num['p']])
    H.evaluate(0,0)
    num=self.num
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(H.output()[0],(q0*tend**6*exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*exp(tend**3/(3*p)))/(3*p**3),9,"Evaluation output mismatch")
    

  def test_issue87(self):
    return # see issue 87
    self.message('IPOPT integration: hessian to p: fwd-over-adjoint on integrator')
    num=self.num
    J=self.qe.jacobian(2)
    J.init()
    J.input(0).set([num['tend']])
    J.input(1).set([num['q0']])
    J.input(2).set([num['p']])
    J.adjSeed(0).set([1])
    J.fwdSeed(0).set([1])
    J.fwdSeed(1).set([1])
    J.fwdSeed(2).set([1])
    # Evaluate
    J.evaluate(1,1)
      
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    print (q0*tend**6*exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*exp(tend**3/(3*p)))/(3*p**3)
    print J.adjSens()
    print J.fwdSens()
    self.assertAlmostEqual(J.adjSens(2)[0],(q0*tend**6*exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*exp(tend**3/(3*p)))/(3*p**3),9,"Evaluation output mismatch")
    
if __name__ == '__main__':
    unittest.main()

