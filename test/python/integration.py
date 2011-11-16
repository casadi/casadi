from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *

scipy_available = True
try:
	import scipy.special
	from scipy.linalg import expm
except:
	scipy_available = False

class Integrationtests(casadiTestCase):

  def setUp(self):
    t=ssym("t")
    q=ssym("q")
    p=ssym("p")
    f = DAE_NUM_IN * [[]]
    f[DAE_T] = t
    f[DAE_Y] = q
    f[DAE_P] = p
    f=SXFunction(f,[q/p*t**2])
    f.init()
    integrator = CVodesIntegrator(f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    integrator.setOption("verbose",True)
    integrator.setOption("t0",0)
    integrator.setOption("tf",2.3)
    integrator.init()
    q0   = MX("q0")
    par  = MX("p")
    qend=integrator.call([q0,par,MX()])[0]
    qe=MXFunction([q0,par],[qend])
    qe.init()
    self.integrator = integrator
    self.qe=qe
    self.qend=qend
    self.q0=q0
    self.par=par
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
    
    f = MXFunction([q0],[qe.call([q0,MX(num['p'])])[0]])
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
    f=SXFunction({'NUM': DAE_NUM_IN, DAE_T: t, DAE_Y: [x,y]},[[z,z]])
    f.init()
    # Pass inputs
    f.setInput(1.0,DAE_T)
    f.setInput([1.0,0.0],DAE_Y)
    # Pass adjoint seeds
    f.setAdjSeed([1.0,0.0])
    # Evaluate with adjoint mode AD
    f.evaluate(0,1)
    # print result
    print f.output()
    print f.adjSens(DAE_Y)
  
  def test_issue92b(self):
    self.message("regression check for issue 92")
    t=SX("t")
    x=SX("x")
    y=SX("y")
    f=SXFunction({'NUM': DAE_NUM_IN, DAE_T: t, DAE_Y: [x,y]},[[x,(1+1e-9)*x]])
    integrator = CVodesIntegrator(f)
    integrator.setOption("t0",0)
    integrator.setOption("tf",1)
    integrator.init()
    # Pass inputs
    integrator.setInput([1,0],INTEGRATOR_X0)
    # Pass adjoint seeds
    integrator.setAdjSeed([1.0,0.0],INTEGRATOR_XF)
    ## Integrate and calculate sensitivities
    integrator.evaluate(0,1)
    # print result
    print integrator.output(INTEGRATOR_XF)
    print integrator.adjSens(INTEGRATOR_X0)
    
  def test_issue92(self):
    self.message("regression check for issue 92")
    t=SX("t")
    x=SX("x")
    var = MX("var",2,1)

    q = [x,SX("problem")]

    dq=[x,x]
    f=SXFunction({'NUM': DAE_NUM_IN, DAE_T: t, DAE_Y: q},[dq])
    f.init()

    integrator = CVodesIntegrator(f)
    integrator.setOption("reltol",1e-12)
    integrator.setOption("t0",0)
    integrator.setOption("tf",1)
    integrator.init()

    qend = integrator.call([var,MX(),MX()])[0]

    f = MXFunction([var],[qend[0]])
    f.init()

    J=Jacobian(f,0)
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
    f=SXFunction([t,q,p,[]],[vertcat([dh ,q[0],dh])])
    f.init()
    
    integrator = CVodesIntegrator(f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    integrator.setOption("verbose",True)
    integrator.setOption("steps_per_checkpoint",10000)
    integrator.setOption("t0",0)
    integrator.setOption("tf",te)

    integrator.init()

    q0   = MX("q0",3,1)
    par  = MX("p",1,1)
    qend=integrator.call([q0,par,MX()])[0]
    qe=MXFunction([q0,par],[qend])
    qe.init()

    #J=self.qe.jacobian(2)
    J=Jacobian(qe,0)
    J.init()
    J.input(0).set(A)
    J.input(1).set(p0)
    J.evaluate()
    outA=J.output().toArray()
    f=SXFunction([t,q,p,[]],[vertcat([dh ,q[0],(1+1e-9)*dh])])
    f.init()
    
    integrator = CVodesIntegrator(f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    integrator.setOption("verbose",True)
    integrator.setOption("steps_per_checkpoint",10000)
    integrator.setOption("t0",0)
    integrator.setOption("tf",te)

    integrator.init()

    q0   = MX("q0",3,1)
    par  = MX("p",1,1)
    qend=integrator.call([q0,par,MX()])[0]
    qe=MXFunction([q0,par],[qend])
    qe.init()

    #J=self.qe.jacobian(2)
    J=Jacobian(qe,0)
    J.init()
    J.input(0).set(A)
    J.input(1).set(p0)
    J.evaluate()
    outB=J.output().toArray()
    print outA-outB
    
  def test_hess(self):
    self.message('CVodes integration: hessian to p: fwd-over-adjoint on integrator')
    num=self.num
    J=self.integrator.jacobian(INTEGRATOR_P,INTEGRATOR_XF)
    J.setOption("number_of_fwd_dir",0)
    J.setOption("number_of_adj_dir",1)
    J.init()
    J.input(INTEGRATOR_X0).set([num['q0']])
    J.input(INTEGRATOR_P).set([num['p']])
    J.adjSeed(INTEGRATOR_XF).set([1])
    # Evaluate
    J.evaluate(0,1)
      
    tend=num['tend']
    q0=num['q0']
    p=num['p']

    self.assertAlmostEqual(J.adjSens(INTEGRATOR_P)[0],(q0*tend**6*exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*exp(tend**3/(3*p)))/(3*p**3),9,"Evaluation output mismatch")

  def test_hess3(self):
    self.message('CVodes integration: hessian to p: Jacobian of integrator.jacobian')
    num=self.num
    J=self.integrator.jacobian(INTEGRATOR_P,INTEGRATOR_XF)
    J.init()
    H=Jacobian(J,INTEGRATOR_P)
    H.setOption("ad_mode","adjoint")
    H.init()
    H.input(INTEGRATOR_X0).set([num['q0']])
    H.input(INTEGRATOR_P).set([num['p']])
    H.evaluate(0,0)
    num=self.num
    tend=num['tend']
    q0=num['q0']
    p=num['p']
    self.assertAlmostEqual(H.output()[0],(q0*tend**6*exp(tend**3/(3*p)))/(9*p**4)+(2*q0*tend**3*exp(tend**3/(3*p)))/(3*p**3),9,"Evaluation output mismatch")

  def test_hess4(self):
    self.message('CVodes integration: hessian to p: Jacobian of integrator.jacobian indirect')
    num=self.num
    J=self.integrator.jacobian(INTEGRATOR_P,INTEGRATOR_XF)
    J.init()
    
    q0=MX("q0")
    p=MX("p")
    dq0=MX("dq0")
    Ji = MXFunction([q0,p,dq0],[J.call([q0,p,dq0])[0]])
    Ji.init()
    H=Jacobian(Ji,1)
    H.setOption("ad_mode","adjoint")
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
    dq0=MX("dq0")
    qe = MXFunction([q0,p,dq0],[self.integrator.call([q0,p,dq0])[0]])
    qe.init()

    JT = MXFunction([q0,p,dq0],[qe.jac(1)[0].T])
    JT.init()
    JT.input(0).set([num['q0']])
    JT.input(1).set([num['p']])
    JT.evaluate(1,0)
    print JT.output()

    H  = Jacobian(JT,1)
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
    dq0=MX("dq0")
    qe = MXFunction([q0,p,dq0],[self.integrator.call([q0,p,dq0])[0]])
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
    f_in = SXMatrixVector(DAE_NUM_IN)
    f_in[DAE_T] = t
    f_in[DAE_Y] = q
    f_in[DAE_P] = p
    f_out = [mul(c.reshape(p,3,3),q)]
    f=SXFunction(f_in,f_out)
    f.init()
    integrator = CVodesIntegrator(f)
    integrator.setOption("steps_per_checkpoint",1000)
    integrator.setOption("t0",0)
    integrator.setOption("tf",te)
    integrator.init()
    q0   = MX("q0",3,1)
    par  = MX("p",9,1)
    qend=integrator.call([q0,par,MX(3,1)])[0]
    qe=integrator.jacobian(INTEGRATOR_P,INTEGRATOR_XF)
    qe.init()
    qe=qe.call([q0,par,MX(3,1)])[0]

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
    f = DAE_NUM_IN * [[]]
    f[DAE_T] = t
    f[DAE_Y] = q
    f[DAE_P] = p
    f=SXFunction(f,[mul(c.reshape(p,3,3),q)])
    f.init()

    integrator = CVodesIntegrator(f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    integrator.setOption("verbose",True)
    integrator.setOption("steps_per_checkpoint",10000)
    integrator.setOption("t0",0)
    integrator.setOption("tf",te)

    integrator.init()

    q0   = MX("q0",3,1)
    par  = MX("p",9,1)
    qend=integrator.call([q0,par,MX(3,1)])[0]
    qe=MXFunction([q0,par],[qend])
    qe.init()
    qendJ=integrator.jacobian(INTEGRATOR_X0,INTEGRATOR_XF)
    qendJ.init()
    qendJ=qendJ.call([q0,par,MX(3,1)])[0]

    qeJ=MXFunction([q0,par],[qendJ])
    qeJ.init()

    qendJ2=integrator.jacobian(INTEGRATOR_X0,INTEGRATOR_XF)
    qendJ2.init()
    qendJ2=qendJ2.call([q0,par,MX(3,1)])[0]

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
    H=Jacobian(qeJ,0,0)
    H.setOption("ad_mode","adjoint")
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

    f=SXFunction([t,q,p,[]],[vertcat([q[1],(p[0]-2*p[1]*cos(2*p[2]))*q[0]])])
    f.init()
    
    integrator = CVodesIntegrator(f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    integrator.setOption("verbose",True)
    integrator.setOption("steps_per_checkpoint",10000)
    integrator.setOption("t0",0)
    integrator.setOption("tf",te)

    integrator.init()

    q0   = MX("q0",2,1)
    par  = MX("p",3,1)
    qend=integrator.call([q0,par,MX(2,1)])[0]
    qe=MXFunction([q0,par],[qend])
    qe.init()
    qendJ=integrator.jacobian(INTEGRATOR_X0,INTEGRATOR_XF)
    qendJ.init()
    qendJ=qendJ.call([q0,par,MX(2,1)])[0]
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
    f=SXFunction([t,q,p,[]],[vertcat([q[1],p[0]+q[1]**2 ])])
    f.init()
    
    integrator = CVodesIntegrator(f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    integrator.setOption("verbose",True)
    integrator.setOption("steps_per_checkpoint",10000)
    integrator.setOption("t0",0)
    integrator.setOption("tf",te)

    integrator.init()

    t0   = MX(0)
    tend = MX(te)
    q0   = MX("q0",2,1)
    par  = MX("p",1,1)
    qend=integrator.call([q0,par,MX(2,1)])[0]
    qe=MXFunction([q0,par],[qend])
    qe.init()
    qendJ=integrator.jacobian(INTEGRATOR_X0,INTEGRATOR_XF)
    qendJ.init()
    qendJ=qendJ.call([q0,par,MX(2,1)])[0]
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
    
    
    Jf=Jacobian(qe,0,0)
    Jf.setOption("ad_mode","adjoint")
    Jf.init()
    Jf.input(0).set(A)
    Jf.input(1).set(p0)
    Jf.evaluate()
    self.checkarray(Jf.output(),Jr,"Jacobian of Nonlin ODE")
    
    
    Jf=Jacobian(qe,0,0)
    Jf.setOption("ad_mode","forward")
    Jf.init()
    Jf.input(0).set(A)
    Jf.input(1).set(p0)
    Jf.evaluate()
    self.checkarray(Jf.output(),Jr,"Jacobian of Nonlin ODE")
    
    # Joel: This is no longer supported: might be a good idea to support again, though
    #qeJ=integrator.jac(INTEGRATOR_X0,INTEGRATOR_XF)
    #qeJ.init()
    #qeJ.input(INTEGRATOR_X0).set(list(A)+[0,1,0,0])
    #qeJ.adjSeed(INTEGRATOR_XF).set([0,0]+[0,1,0,0])
    #qeJ.evaluate(0,1)
    #print qeJ.output()
    #print qeJ.adjSens(INTEGRATOR_X0)
    
    Jr = matrix([[(sqrt(p0)*(te*yc0**2-yc0+p0*te)*tan(arctan(yc0/sqrt(p0))+sqrt(p0)*te)+yc0**2)/(2*p0*yc0**2+2*p0**2)],[(sqrt(p0)*((te*yc0**2-yc0+p0*te)*tan(arctan(yc0/sqrt(p0))+sqrt(p0)*te)**2+te*yc0**2-yc0+p0*te)+(yc0**2+p0)*tan(arctan(yc0/sqrt(p0))+sqrt(p0)*te))/(sqrt(p0)*(2*yc0**2+2*p0))]])  
    
    Jf=Jacobian(qe,1,0)
    Jf.setOption("ad_mode","adjoint")
    Jf.init()
    Jf.input(0).set(A)
    Jf.input(1).set(p0)
    Jf.evaluate()
    self.checkarray(Jf.output(),Jr,"Jacobian of Nonlin ODE")
    
    Jf=Jacobian(qe,1,0)
    Jf.setOption("ad_mode","forward")
    Jf.init()
    Jf.input(0).set(A)
    Jf.input(1).set(p0)
    Jf.evaluate()
    self.checkarray(Jf.output(),Jr,"Jacobian of Nonlin ODE")
    
    qendJ=integrator.jacobian(INTEGRATOR_P,INTEGRATOR_XF)
    qendJ.init()
    qendJ=qendJ.call([q0,par,MX(2,1)])[0]
    qeJ=MXFunction([q0,par],[qendJ])
    qeJ.init()

    qeJ.input(0).set(A)
    qeJ.input(1).set(p0)
    qeJ.evaluate()
    
    self.checkarray(qeJ.output(),Jr,"jacobian of Nonlin ODE")
    
    
    
    
    qeJf=MXFunction([q0,par],[vec(qeJ.call([q0,par])[0])])
    qeJf.init()
    
    H=Jacobian(qeJf,0,0)
    H.setOption("ad_mode","adjoint")
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
    #qeJ=integrator.jac(INTEGRATOR_X0,INTEGRATOR_XF)
    #qeJ.init()
    #qeJ.input(INTEGRATOR_X0).set(list(A)+[0,1,0,0])
    #qeJ.adjSeed(INTEGRATOR_XF).set([0,0]+[0,1,0,0])
    #qeJ.evaluate(0,1)
    #print qeJ.output()
    #print qeJ.adjSens(INTEGRATOR_X0)
    
  def test_hessian2D(self):
    self.message("hessian")
    N=2

    x0_ = DMatrix([1,0.1])
    A_  = DMatrix([[3,1],[0.74,4]])

    A = ssym("A",N,N)
    x = ssym("x",N)

    ode = SXFunction({'NUM':DAE_NUM_IN, DAE_Y: x, DAE_P: A},[mul(A,x)])
    I = CVodesIntegrator(ode)
    I.init()
    I.setOption('reltol',1e-12)
    I.input(INTEGRATOR_X0).set(x0_)
    I.input(INTEGRATOR_P).set(A_)
    I.evaluate()

    q0=MX("q0",N)
    p=MX("p",N*N)
    dq0=MX("dq0",N)
    qe = MXFunction([q0,p,dq0],[I.call([q0,p,dq0])[0]])
    qe.init()

    JT = MXFunction([q0,p,dq0],[qe.jac(1)[0].T])
    JT.init()

    H  = Jacobian(JT,1)
    H.init()
    H.input(0).set(x0_)
    H.input(1).set(vec(A_))
    H.evaluate()

    H1 = DMatrix(H.output())
    
    H = qe.hessian(1)
    H.init()
    H.input(0).set(x0_)
    H.input(1).set(vec(A_))
    H.evaluate()
    H2 = DMatrix(H.output())
    
    self.checkarray(H1,H2,"hessian")
    print H.output()

    
if __name__ == '__main__':
    unittest.main()

