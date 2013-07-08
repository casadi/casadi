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
import unittest
from types import *
from helpers import *
from copy import deepcopy

scipy_available = True
try:
	from scipy.sparse import csr_matrix
except:
	scipy_available = False

def checkarray(self,zr,zt,name):
    if len(zr.shape)==1 and (zt.shape[0]==1 or zt.shape[1]==1) and zr.shape[0]==zt.shape[1]*zt.shape[0]:
      zr=reshape(zr,(zt.shape));
    self.assertEqual(zt.shape[0],zr.shape[0],"%s dimension error. Got %s, expected %s" % (name,str(zt.shape),str(zr.shape)))
    self.assertEqual(len(zt.shape),len(zr.shape),"%s dimension error. Got %s, expected %s" % (name,str(zt.shape),str(zr.shape)))
    self.assertEqual(zt.shape[1],zr.shape[1],"%s dimension error. Got %s, expected %s" % (name,str(zt.shape),str(zr.shape)))
    for i in range(zr.shape[0]):
      for j in range(zr.shape[1]):
        self.assertAlmostEqual(zt[i,j],zr[i,j],10,"%s evaluation error. %s <-> %s" % (name, str(zt),str(zr)))

def checkMXoperations(self,ztf,zrf,name):
    x = MX("x",1,3)
    z=vertcat([x*(i+1) for i in range(8)])
    f = MXFunction([x],[ztf(z)])
    f.init()
    L=[1,2,3]
    f.setInput(L,0)
    f.evaluate()
    zt = f.output(0).toArray()
    zr = array([[L[0]*(i+1),L[1]*(i+1),L[2]*(i+1)] for i in range(8)])
    checkarray(self,zrf(zr),zt,name)
    return (zt,zrf(zr))

def checkMXoperations2(self,ztf,zrf,name):
    x = MX("x",3,1)
    z = horzcat([x*i for i in range(8)])
    f = MXFunction([x],[ztf(z)])
    f.init()
    L=[1,2,3]
    f.setInput(L,0)
    f.evaluate()
    zt = f.output(0).toArray()
    zr = array([[L[0]*i,L[1]*i,L[2]*i] for i in range(8)]).T
    checkarray(self,zrf(zr),zt,name)
    return zt

def checkMXoperations3(self,ztf,zrf,name):
    x = MX("x",3,1)
    p = horzcat([x[0,0],x[1,0],x[2,0]])
    z = vertcat([p*i for i in range(8)])
    f = MXFunction([x],[ztf(z)])
    f.init()
    L=[1,2,3]
    f.setInput(L,0)
    f.evaluate()
    zt = f.output(0).toArray()
    zr = array([[L[0]*i,L[1]*i,L[2]*i] for i in range(8)])
    checkarray(self,zrf(zr),zt,name)
    return (zt,zrf(zr))
    
class MXtests(casadiTestCase):

  def setUp(self):
    self.pool=FunctionPool()
    self.pool.append(lambda x: sqrt(x[0]),sqrt,"sqrt")
    self.pool.append(lambda x: sin(x[0]),sin,"sin")
    self.pool.append(lambda x: cos(x[0]),cos,"cos")
    self.pool.append(lambda x: tan(x[0]),tan,"tan")
    self.pool.append(lambda x: arctan(x[0]),arctan,"arctan")
    self.pool.append(lambda x: arcsin(x[0]),arcsin,"arcsin")
    self.pool.append(lambda x: arccos(x[0]),arccos,"arccos")
    self.pool.append(lambda x: exp(x[0]),exp,"exp")
    self.pool.append(lambda x: log(x[0]),log,"log")
    self.pool.append(lambda x: x[0]**0,lambda x : x**0,"x^0",flags={'nozero'})
    self.pool.append(lambda x: x[0]**1,lambda x : x**1,"^1")
    self.pool.append(lambda x: x[0]**(-2),lambda x : x**(-2),"^-2",flags={'nozero'})
    self.pool.append(lambda x: x[0]**(0.3),lambda x : x**(0.3),"^0.3")
    self.pool.append(lambda x: floor(x[0]),floor,"floor")
    self.pool.append(lambda x: ceil(x[0]),ceil,"ceil")
    self.Jpool=FunctionPool()
    self.Jpool.append(lambda x: sqrt(x[0]),lambda x:diag(1/(2.0*sqrt(x))),"sqrt")
    self.Jpool.append(lambda x: sin(x[0]),lambda x:diag(cos(x)),"sin")
    self.Jpool.append(lambda x: cos(x[0]),lambda x:diag(-sin(x)),"cos")
    self.Jpool.append(lambda x: tan(x[0]),lambda x:diag(1.0/cos(x)**2),"tan")
    self.Jpool.append(lambda x: arctan(x[0]),lambda x:diag( 1.0/(x**2+1)),"arctan")
    self.Jpool.append(lambda x: arcsin(x[0]),lambda x:diag( 1.0/sqrt(1-x**2)),"arcsin")
    self.Jpool.append(lambda x: arccos(x[0]),lambda x: diag(-1.0/sqrt(1-x**2)),"arccos")
    self.Jpool.append(lambda x: exp(x[0]),lambda x: diag(exp(x)),"exp")
    self.Jpool.append(lambda x: log(x[0]),lambda x: diag(1.0/x),"log")
    self.Jpool.append(lambda x: x[0]**0,lambda x :diag(zeros(x.shape)),"x^0")
    self.Jpool.append(lambda x: x[0]**1,lambda x : diag(ones(x.shape)),"^1")
    self.Jpool.append(lambda x: x[0]**(-2),lambda x : diag(-2.0/x**3),"^-2")
    self.Jpool.append(lambda x: x[0]**(0.3),lambda x :diag( 0.3/x**0.7),"^0.3")
    self.matrixpool=FunctionPool()
    #self.matrixpool.append(lambda x: norm_2(x[0]),linalg.norm,"norm_2")
    #self.matrixpool.append(lambda x: norm_1(x[0]),lambda x: sum(sum(abs(x))),"norm_1")
    #self.matrixpool.append(lambda x: norm_inf(x[0]),lambda x: abs(matrix(x)).max(),"norm_inf")
    self.matrixbinarypool=FunctionPool()
    self.matrixbinarypool.append(lambda a: a[0]+a[1],lambda a: a[0]+a[1],"Matrix+Matrix")
    self.matrixbinarypool.append(lambda a: a[0]-a[1],lambda a: a[0]-a[1],"Matrix-Matrix")
    self.matrixbinarypool.append(lambda a: a[0]*a[1],lambda a: a[0]*a[1],"Matrix*Matrix")
    #self.matrixbinarypool.append(lambda a: inner_mul(a[0],trans(a[1])),lambda a: dot(a[0].T,a[1]),name="inner_mul(Matrix,Matrix)") 
    self.matrixbinarypool.append(lambda a: mul(a[0],trans(a[1])),lambda a: dot(a[0],a[1].T),"mul(Matrix,Matrix.T)")

  def test_indirection(self):
    self.message("MXFunction indirection")
    x=MX("x",2,1)
    y=MX("y",2,1)
    z=MX("z",2,1)
    
    xn = array([2.3,1.3])
    yn = array([7.3,4.6])
    zn = array([12,7.4])
    
    f=MXFunction([x,y,z],[x+2*y+3*z])
    f.init()
    self.message(":simple indirection")
    g=MXFunction([x,y,z],f.call([x,y,z]))
    g.init()
    g.setInput(xn,0)
    g.setInput(yn,1)
    g.setInput(zn,2)
    g.setFwdSeed(xn,0) # okay, I'm just lazy comming up with more numbers
    g.setFwdSeed(yn,1)
    g.setFwdSeed(zn,2)
    g.evaluate(1,0)
    self.checkarray(g.getOutput(),xn+2*yn+3*zn,"MXFunction indirection");
    self.checkarray(g.getFwdSens(),array([52.9,32.7]),"MXFunction indirection");
    
    g=MXFunction([x,y,z],f.call([vertcat([x[0],x[1]]),y,z]))
    g.init()
    g.setInput(xn,0)
    g.setInput(yn,1)
    g.setInput(zn,2)
    g.setFwdSeed(xn,0)
    g.setFwdSeed(yn,1)
    g.setFwdSeed(zn,2)
    g.evaluate(1,0)

    self.checkarray(g.getOutput(),xn+2*yn+3*zn,"MXFunction indirection");
    self.checkarray(g.getFwdSens(),array([52.9,32.7]),"MXFunction indirection");
    
    self.message(":double output flipover")
    h=MXFunction([x,y,z],f.call([vertcat([y[0],x[1]]),vertcat([x[0],y[1]]),z]))
    h.init()
    h=MXFunction([x,y,z],h.call([vertcat([y[0],x[1]]),vertcat([x[0],y[1]]),z]))
    # h should be identical to g now
    h.init()
    h.setInput(xn,0)
    h.setInput(yn,1)
    h.setInput(zn,2)
    h.setFwdSeed(xn,0)
    h.setFwdSeed(yn,1)
    h.setFwdSeed(zn,2)
    h.evaluate(1,0)
    self.checkarray(h.getOutput(),xn+2*yn+3*zn,"MXFunction indirection");
    self.checkarray(h.getFwdSens(),array([52.9,32.7]),"MXFunction indirection");
    
    self.message(":double input flipover")
    h=MXFunction([x,y,z],f.call([y,x,z]))
    h.init()
    h=MXFunction([x,y,z],h.call([y,x,z]))
    h.init()
    h.setInput(xn,0)
    h.setInput(yn,1)
    h.setInput(zn,2)
    h.setFwdSeed(xn,0)
    h.setFwdSeed(yn,1)
    h.setFwdSeed(zn,2)
    h.evaluate(1,0)
    self.checkarray(h.getOutput(),xn+2*yn+3*zn,"MXFunction indirection");
    self.checkarray(h.getFwdSens(),array([52.9,32.7]),"MXFunction indirection");
    
    return # uncomplete calls are not supported
    self.message(":uncomplete call")
    h=MXFunction([x,z],f.call([x,y,z]))
    h.init()
    h=MXFunction([x,y,z],h.call([x,z]))
    h.init()
    h.setInput(xn,0)
    h.setInput(yn,1)
    h.setInput(zn,2)
    h.setFwdSeed(xn,0)
    h.setFwdSeed(yn,1)
    h.setFwdSeed(zn,2)
    h.evaluate(1,0)
    self.checkarray(h.getOutput(),xn+2*yn+3*zn,"MXFunction indirection");
    self.checkarray(h.getFwdSens(),array([52.9,32.7]),"MXFunction indirection");
    
  def test_MX1(self):
    self.message("MX constructor")
    x = MX("x",2,3)
    self.assertEqual(x.size1(),2,"MX fails to indicate its size1")
    self.assertEqual(x.size2(),3,"MX fails to indicate its size2")

  def test_MXvertcat(self):
    self.message("MX vertcat")
    x = MX("x",1,3)
    y = MX("y",1,3)
    z=vertcat((x,y))
    self.assertEqual(z.size1(),2,"MX fails to indicate its size1")
    self.assertEqual(z.size2(),3,"MX fails to indicate its size2")

  def test_MXFunction1(self):
    self.message("MXFunction single input, single output")
    # check if x->2*x
    # evaluates correctly for x=3
    x = MX("x")
    y = 2*x
    f = MXFunction([x],[y])
    self.assertEqual(f.getNumInputs(),1,"MXFunction fails to indicate correct number of inputs")
    self.assertEqual(f.getNumOutputs(),1,"MXFunction fails to indicate correct number of outputs")
    f.init()
    f.setInput(3,0);
    f.evaluate()
    yt = tuple(f.output().data())
    self.assertEqual(type(yt),TupleType,"Output of MXFunction is expected to be tuple of floats")
    self.assertEqual(len(yt),1,"Output of MXFunction was tuple of floats, as expected, but length is incorrect.")
    y=yt[0]
    self.assertEqual(type(y),float,"Output of MXFunction is expected to be tuple of floats")
    self.assertAlmostEqual(y, 2*3,10)

  def test_MXfunction2(self):
    self.message("MXFunction multi input, multi output")
      # check if [x,y]->[y+x,y*x]
    # evaluates correctly for x=3,y=7
    x = MX("x")
    y = MX("y")
    f = MXFunction([x,y],[x+y,y*x])
    self.assertEqual(f.getNumInputs(),2,"MXFunction fails to indicate correct number of inputs")
    self.assertEqual(f.getNumOutputs(),2,"MXFunction fails to indicate correct number of outputs")

    f.init()
    f.setInput(3,0);
    f.setInput(7,1);
    f.evaluate()
    zt1 = tuple(f.output(0).data())
    zt2 = tuple(f.output(1).data())
    self.assertEqual(type(zt1),TupleType,"Output of MXFunction is expected to be tuple of floats")
    self.assertEqual(type(zt2),TupleType,"Output of MXFunction is expected to be tuple of floats")
    self.assertEqual(len(zt1),1,"Output of MXFunction was tuple of floats, as expected, but length is incorrect.")
    self.assertEqual(len(zt2),1,"Output of MXFunction was tuple of floats, as expected, but length is incorrect.")
    z1=zt1[0]
    z2=zt2[0]
    self.assertEqual(type(z1),float,"Output of MXFunction is expected to be tuple of floats")
    self.assertEqual(type(z2),float,"Output of MXFunction is expected to be tuple of floats")
    self.assertAlmostEqual(z2, 21,10)
    self.assertAlmostEqual(z1, 10,10)      



  def test_MXfunction3(self):
    self.message("MXFunction single input, multi output (1)")
    # check if [x,y]->[y+x,y*x]
    # evaluates correctly for x=3,y=7
    # now with single input, multi output
    xy = MX("xy",2)
    f = MXFunction([xy],[xy[0]+xy[1],xy[0]*xy[1]])
    self.assertEqual(f.getNumInputs(),1,"MXFunction fails to indicate correct number of inputs")
    self.assertEqual(f.getNumOutputs(),2,"MXFunction fails to indicate correct number of outputs")
    f.init()
    f.setInput([3,7],0);
    f.evaluate()
    zt1 = tuple(f.output(0).data())
    zt2 = tuple(f.output(1).data())
    self.assertEqual(type(zt1),TupleType,"Output of MXFunction is expected to be tuple of floats")
    self.assertEqual(type(zt2),TupleType,"Output of MXFunction is expected to be tuple of floats")
    self.assertEqual(len(zt1),1,"Output of MXFunction was tuple of floats, as expected, but length is incorrect.")
    self.assertEqual(len(zt2),1,"Output of MXFunction was tuple of floats, as expected, but length is incorrect.")
    z1=zt1[0]
    z2=zt2[0]
    self.assertEqual(type(z1),float,"Output of MXFunction is expected to be tuple of floats")
    self.assertEqual(type(z2),float,"Output of MXFunction is expected to be tuple of floats")
    self.assertAlmostEqual(z2, 21,10)
    self.assertAlmostEqual(z1, 10,10)

  def test_MXfunction3b(self):
    self.message("MXFunction single input, multi output (2)")
    # check if [x,y]->[y+x,y*x]
    # evaluates correctly for x=3,y=7
    # now with single input, multi output
    xy = MX("xy",1,2)
    f = MXFunction([xy],[xy[0]+xy[1],xy[0]*xy[1]])
    self.assertEqual(f.getNumInputs(),1,"MXFunction fails to indicate correct number of inputs")
    self.assertEqual(f.getNumOutputs(),2,"MXFunction fails to indicate correct number of outputs")
    f.init()
    f.setInput([3,7],0);
    f.evaluate()
    zt1 = f.output(0).toArray()
    zt2 = f.output(1).toArray()
    
    self.assertEqual(type(zt1),ndarray,"Output of MXFunction is expected to be numpy.ndarray")
    self.assertEqual(zt1.shape[0],1,"Output of MXFunction is of wrong shape.")
    self.assertEqual(zt1.shape[1],1,"Output of MXFunction is of wrong shape.")
    
    self.assertEqual(type(zt2),ndarray,"Output of MXFunction is expected to be numpy.ndarray")
    self.assertEqual(zt2.shape[0],1,"Output of MXFunction is of wrong shape.")
    self.assertEqual(zt2.shape[1],1,"Output of MXFunction is of wrong shape.")
    
    z1=zt1[0,0]
    z2=zt2[0,0]
    self.assertEqual(type(z1),float64,"Output of MXFunction is expected to be numpy.ndarray of floats")
    self.assertEqual(type(z2),float64,"Output of MXFunction is expected to be numpy.ndarray of floats")
    self.assertAlmostEqual(z2, 21,10)
    self.assertAlmostEqual(z1, 10,10)

  def test_MXfunction4(self):
    self.message("MXFunction single input, single output , using vertcat")
    # check if [x,y]->[y+x,y*x]
    # evaluates correctly for x=3,y=7
    # now with single input, single output
    xy = MX("xy",2)
    z=vertcat([xy[0]+xy[1],xy[0]*xy[1]])
    f = MXFunction([xy],[z])
    self.assertEqual(f.getNumInputs(),1,"MXFunction fails to indicate correct number of inputs")
    self.assertEqual(f.getNumOutputs(),1,"MXFunction fails to indicate correct number of outputs")
    f.init()
    f.setInput([3,7],0);
    f.evaluate()
    zt=f.output(0).toArray()
    self.assertEqual(type(zt),ndarray,"Output of MXFunction is expected to be numpy.ndarray")
    self.assertEqual(zt.shape[0],2,"Output of MXFunction is of wrong shape.")
    self.assertEqual(zt.shape[1],1,"Output of MXFunction is of wrong shape.")
    z1=zt[0,0]
    z2=zt[1,0]
    self.assertEqual(type(z1),float64,"Output of MXFunction is expected to be numpy.ndarray of floats")
    self.assertEqual(type(z2),float64,"Output of MXFunction is expected to be numpy.ndarray of floats")
    self.assertAlmostEqual(z2, 21,10)
    self.assertAlmostEqual(z1, 10,10)

  def test_MXfunction5(self):
    self.message("MXFunction single input, single output , using horzcat")
    # check if [x,y]->[y+x,y*x]
    # evaluates correctly for x=3,y=7
    # now with single input, single output
    xy = MX("xy",2)
    z=horzcat([xy[0]+xy[1],xy[0]*xy[1]])
    f = MXFunction([xy],[z])
    self.assertEqual(f.getNumInputs(),1,"MXFunction fails to indicate correct number of inputs")
    self.assertEqual(f.getNumOutputs(),1,"MXFunction fails to indicate correct number of outputs")
    f.init()
    f.setInput([3,7],0);
    f.evaluate()
    zt = f.output(0).toArray()
    self.assertEqual(type(zt),ndarray,"Output of MXFunction is expected to be numpy.ndarray")
    self.assertEqual(zt.shape[0],1,"Output of MXFunction is of wrong shape.")
    self.assertEqual(zt.shape[1],2,"Output of MXFunction is of wrong shape.")
    z1=zt[0,0]
    z2=zt[0,1]
    self.assertEqual(type(z1),float64,"Output of MXFunction is expected to be numpy.ndarray of floats")
    self.assertEqual(type(z2),float64,"Output of MXFunction is expected to be numpy.ndarray of floats")
    self.assertAlmostEqual(z2, 21,10)
    self.assertAlmostEqual(z1, 10,10)
    
  def test_issue83(self):
    x=MX("x")
    y=MX("y")

    z = x + y

    f = MXFunction([x,y],[z])
    f.init()

    [fc] = f.call([MX(3),y])

    g = MXFunction([y],[fc])
    g.init()
    g.setInput([7])
    g.evaluate()

    self.assertAlmostEqual(g.getOutput()[0],10,10,"issue #83")
    
    [fc] = f.call([x,MX(7)])

    g = MXFunction([x],[fc])
    g.init()
    g.setInput([3])
    g.evaluate()

    self.assertAlmostEqual(g.getOutput()[0],10,10,"issue #83")
        
  def test_identitySX(self):
    self.message("identity SXFunction")
    x = SX("x")
    f = SXFunction([x],[x])
    f.init()
    f.setInput([3],0)
    f.evaluate()
    self.assertAlmostEqual(f.getOutput(0)[0,0], 3,10)

  def test_identityMX(self):
    self.message("identity MXFunction")
    x = MX("x")
    f  = MXFunction([x],[x])
    f.init()
    f.setInput([3],0)
    f.evaluate()
    self.assertAlmostEqual(f.getOutput(0)[0,0], 3,10)
    
  def test_MXorder(self):
    self.message("MXFunction order of non-zero elements")
    x = MX("x",2,3)
    f = MXFunction([x],[x+x])

    self.assertEqual(f.getNumInputs(),1,"MXFunction fails to indicate correct number of inputs")
    self.assertEqual(f.getNumOutputs(),1,"MXFunction fails to indicate correct number of outputs")
    f.init()
    L=[1,2,3,4,5,6]
    f.setInput(L,0)
    f.evaluate()
    zt = f.output(0).toArray()
    self.assertEqual(zt.shape[0],2,"Output of MXFunction is of wrong shape.")
    self.assertEqual(zt.shape[1],3,"Output of MXFunction is of wrong shape.")
      
    Lr=reshape(L,(2,3))
    for i in range(2):
      for j in range(3):
        self.assertAlmostEqual(Lr[i,j]*2, zt[i,j],10)
    
  def test_trans(self):
    self.message("trans")
    a = MX(0,1)
    b = trans(a)
    self.assertEquals(b.size1(),1)
    self.assertEquals(b.size2(),0)
    
  def test_MXtrans(self):
    self.message("trans(MX)")
    x = MX("x",2,3)
    z=trans(x)
    self.assertEqual(z.size1(),3,"Flatten returns MX of wrong dimension")
    self.assertEqual(z.size2(),2,"Flatten returns MX of wrong dimension")
    f = MXFunction([x],[z])
    self.assertEqual(f.getNumInputs(),1,"MXFunction fails to indicate correct number of inputs")
    self.assertEqual(f.getNumOutputs(),1,"MXFunction fails to indicate correct number of outputs")
    f.init()
    L=[1,2,3,4,5,6]
    f.setInput(L,0)
    f.evaluate()
    zt = f.output(0).toArray()
    
    ztr=reshape(zt,(3,2))
    Lr=reshape(L,(2,3))
    for i in range(2):
      for j in range(3):
        self.assertAlmostEqual(Lr[i,j], ztr[j,i],10)
      
  def test_MXvec(self):

    u = DMatrix([[10*j+i for i in range(3)] for j in range(4) ])

    U = msym("u",u.shape)

    f = MXFunction([U],[vec(U)])
    f.init()
    f.setInput(u)
    f.evaluate()
    
    self.checkarray(vec(u),f.getOutput(),"vec")
    
  def test_MXvecNZ(self):

    u = DMatrix(4,3)
    u[0,0] = 7
    u[1,1] = 8
    u[2,2] = 6
    u[1,0] = 9
    u[0,1] = 11
    
    U = msym("u",u.sparsity())

    f = MXFunction([U],[vecNZ(U)])
    f.init()
    f.setInput(u)
    f.evaluate()
    
    self.checkarray(vecNZ(u),f.getOutput(),"vec")

  def test_MXreshape(self):
    self.message("reshape(MX)")
    x = MX("x",2,3)
    z=c.reshape(x,(1,6))
    self.assertEqual(z.size1(),1,"Flatten returns MX of wrong dimension")
    self.assertEqual(z.size2(),6,"Flatten returns MX of wrong dimension")
    f = MXFunction([x],[z])
    self.assertEqual(f.getNumInputs(),1,"MXFunction fails to indicate correct number of inputs")
    self.assertEqual(f.getNumOutputs(),1,"MXFunction fails to indicate correct number of outputs")
    f.init()
    L=[1,2,3,4,5,6]
    f.setInput(L,0)
    f.evaluate()
    zt = f.output(0).toArray()
    for i in range(len(L)):
      self.assertAlmostEqual(L[i], zt[0,i],10)
  
  def test_MXcompose(self):
    self.message("compositions of flatten, trans, reshape with vertcat")
    checkMXoperations(self,lambda x: x,lambda x: x,'vertcat')
    checkMXoperations(self,lambda x: trans(x),lambda x: x.T,'trans(vertcat)')
    checkMXoperations(self,lambda x: trans(trans(x)),lambda x: x,'trans(trans(vertcat))')
    checkMXoperations(self,lambda x: vec(trans(x)),lambda x: reshape(x,(prod(x.shape),1)),'vec(trans(vertcat))')
    checkMXoperations(self,lambda x: trans(vec(x)),lambda x: reshape(x.T,(prod(x.shape),1)).T,'vec(trans(vertcat))')
    checkMXoperations(self,lambda x: c.reshape(x,(4,6)),lambda x: reshape(x,(4,6)),'reshape(vertcat)')
    checkMXoperations(self,lambda x: c.reshape(trans(x),(4,6)),lambda x: reshape(x.T,(4,6)),'reshape(trans(vertcat))') 
    checkMXoperations(self,lambda x: trans(c.reshape(x,(4,6))),lambda x: reshape(x,(4,6)).T,'trans(reshape(vertcat))') 

  def test_MXcompose2(self):
    self.message("compositions of flatten, trans, reshape with horzcat")
    checkMXoperations2(self,lambda x: x,lambda x: x,'horzcat')
    checkMXoperations2(self,lambda x: trans(x),lambda x: x.T,'trans(horzcat)')
    checkMXoperations2(self,lambda x: trans(trans(x)),lambda x: x,'trans(trans(horzcat))')
    checkMXoperations2(self,lambda x: vec(trans(x)),lambda x: reshape(x,(prod(x.shape),1)),'vec(trans(horzcat))')
    checkMXoperations2(self,lambda x: trans(vec(x)),lambda x: reshape(x.T,(prod(x.shape),1)).T,'vec(trans(horzcat))')
    checkMXoperations2(self,lambda x: c.reshape(x,(4,6)),lambda x: reshape(x,(4,6)),'reshape(horzcat)')
    checkMXoperations2(self,lambda x: c.reshape(trans(x),(4,6)),lambda x: reshape(x.T,(4,6)),'reshape(trans(horzcat))') 
    checkMXoperations2(self,lambda x: trans(c.reshape(x,(4,6))),lambda x: reshape(x,(4,6)).T,'trans(reshape(horzcat))') 

  def test_MXcompose3(self):
    self.message("compositions of flatten, trans, reshape with vertcat")
    checkMXoperations3(self,lambda x: x,lambda x: x,'snippet')
    checkMXoperations3(self,lambda x: trans(x),lambda x: x.T,'trans(snippet)')
    checkMXoperations3(self,lambda x: trans(trans(x)),lambda x: x,'trans(trans(snippet))')
    checkMXoperations3(self,lambda x: vec(trans(x)),lambda x: reshape(x,(prod(x.shape),1)),'vec(trans(snippet))')
    checkMXoperations3(self,lambda x: trans(vec(x)),lambda x: reshape(x.T,(prod(x.shape),1)).T,'vec(trans(snippet))')
    checkMXoperations3(self,lambda x: c.reshape(x,(4,6)),lambda x: reshape(x,(4,6)),'reshape(snippet)')
    checkMXoperations3(self,lambda x: c.reshape(trans(x),(4,6)),lambda x: reshape(x.T,(4,6)),'reshape(trans(snippet))') 
    checkMXoperations3(self,lambda x: trans(c.reshape(x,(4,6))),lambda x: reshape(x,(4,6)).T,'trans(reshape(snippet))') 

  def test_MXcompose4(self):
    self.message("compositions of horzcat + vertcat")
    checkMXoperations(self,lambda x: vertcat([x]),lambda x: x,'vertcat(vertcat)')
    checkMXoperations(self,lambda x: vertcat([x,x*2]),lambda x: vstack((x,x*2)),'vertcat(vertcat,vertcat)')
    checkMXoperations(self,lambda x: horzcat([x]),lambda x: x,'horzcat(vertcat)')
    checkMXoperations(self,lambda x: horzcat([x,x*2]),lambda x: hstack((x,x*2)),'horzcat(vertcat,vertcat)')
    
    checkMXoperations2(self,lambda x: vertcat([x]),lambda x: x,'vertcat(horzcat)')
    checkMXoperations2(self,lambda x: vertcat([x,x*2]),lambda x: vstack((x,x*2)),'vertcat(horzcat,horzcat)')
    checkMXoperations2(self,lambda x: horzcat([x]),lambda x: x,'horzcat(horzcat)')
    checkMXoperations2(self,lambda x: horzcat([x,x*2]),lambda x: hstack((x,x*2)),'horzcat(horzcat,horzcat)')
    
    checkMXoperations3(self,lambda x: vertcat([x]),lambda x: x,'vertcat(snippet)')
    checkMXoperations3(self,lambda x: vertcat([x,x*2]),lambda x: vstack((x,x*2)),'vertcat(snippet,snippet)')
    checkMXoperations3(self,lambda x: horzcat([x]),lambda x: x,'horzcat(snippet)')
    checkMXoperations3(self,lambda x: horzcat([x,x*2]),lambda x: hstack((x,x*2)),'horzcat(snippet,snippet)')
    
  def test_MXslicingnew(self):
    self.message("MX slicing new")

    self.message(":dense")
    x = MX("x",3,2)
    x0=array([[0.738,0.2],[ 0.1,0.39 ],[0.99,0.999999]])
    self.numpyEvaluationCheck(lambda x:x[0][1,0],lambda x: x[1,0],[x],x0,'x[1,0]')
    self.numpyEvaluationCheck(lambda x:x[0][0,1],lambda x: x[0,1],[x],x0,'x[0,1]')
    self.numpyEvaluationCheck(lambda x:x[0][0,-1],lambda x: x[0,-1],[x],x0,'x[0,-1]')
    self.numpyEvaluationCheck(lambda x: x[0][:,0], lambda x: matrix(x)[:,0],[x],x0,name="x[:,0]")
    self.numpyEvaluationCheck(lambda x: x[0][:,1], lambda x: matrix(x)[:,1],[x],x0,name="x[:,1]")
    self.numpyEvaluationCheck(lambda x: x[0][1,:], lambda x: matrix(x)[1,:],[x],x0,name="x[1,:]")
    self.numpyEvaluationCheck(lambda x: x[0][0,:], lambda x: matrix(x)[0,:],[x],x0,name="x[0,:]")
    self.numpyEvaluationCheck(lambda x: x[0][-1,:], lambda x: matrix(x)[-1,:],[x],x0,name="x[-1,:]")
    self.numpyEvaluationCheck(lambda x: x[0][:,-2], lambda x: matrix(x)[:,-2],[x],x0,name="x[:,-2]")
    self.numpyEvaluationCheck(lambda x: x[0][0:-2,0:-1], lambda x: matrix(x)[0:-2,0:-1],[x],x0,name="x[0:-2,0:-1]")
    self.numpyEvaluationCheck(lambda x: x[0][0:2,0:2], lambda x: matrix(x)[0:2,0:2],[x],x0,name="x[0:2,0:2]")
    self.numpyEvaluationCheck(lambda x: x[0][[0,1],0:2], lambda x: matrix(x)[[0,1],0:2],[x],x0,name="x[[0,1],0:2]")
    self.numpyEvaluationCheck(lambda x: x[0][[0,2,3]], lambda x: matrix([x[0,0],x[1,0],x[1,1]]).T,[x],x0,name="x[[0,2,3]]")

    self.numpyEvaluationCheck(lambda x: x[0][1], lambda x: matrix(x.ravel()[1]).T,[x],x0,name="x[1] on dense matrix")
    self.numpyEvaluationCheck(lambda x: x[0][-1], lambda x: matrix(x.ravel()[-1]).T,[x],x0,name="x[-1] on dense matrix")
    self.numpyEvaluationCheck(lambda x: x[0][[0,1],0:1],lambda x: x[[0,1],0:1],[x],x0,name='x[:,0:1]')
    self.numpyEvaluationCheck(lambda x: x[0][0:1,[0,1]],lambda x: x[0:1,[0,1]],[x],x0,name='x[0:1,:]')
    
    self.message(":sparse")
      
    sp=CRSSparsity(3,4,[1,2,1],[0,2,2,3])
    x=MX("X",sp)
    sx0=[0.738,0.39,0.99]
    x0=DMatrix(3,4,[1,2,1],[0,2,2,3],[0.738,0.39,0.99]).toArray()
    self.numpyEvaluationCheck(lambda x: x[0][0,0], lambda x: matrix(x)[0,0],[x],x0,name="x[0,0]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][1,0], lambda x: matrix(x)[1,0],[x],x0,name="x[1,0]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][0,1], lambda x: matrix(x)[0,1],[x],x0,name="x[1,0]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][0,-1], lambda x: matrix(x)[0,-1],[x],x0,name="x[0,-1]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][:,0], lambda x: matrix(x)[:,0],[x],x0,name="x[:,0]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][:,1], lambda x: matrix(x)[:,1],[x],x0,name="x[:,1]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][1,:], lambda x: matrix(x)[1,:],[x],x0,name="x[1,:]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][0,:], lambda x: matrix(x)[0,:],[x],x0,name="x[0,:]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][-1,:], lambda x: matrix(x)[-1,:],[x],x0,name="x[-1,:]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][:,-2], lambda x: matrix(x)[:,-2],[x],x0,name="x[:,-2]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][0:-2,0:-1], lambda x: matrix(x)[0:-2,0:-1],[x],x0,name="x[0:-2,0:-1]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][0:2,0:2], lambda x: matrix(x)[0:2,0:2],[x],x0,name="x[0:2,0:2]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][[0,1],0:2], lambda x: matrix(x)[[0,1],0:2],[x],x0,name="x[[0,1],0:2]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][[2,1]], lambda x: matrix([x[2,1],x[0,2]]).T,[x],x0,name="x[[2,1]]")
    self.numpyEvaluationCheck(lambda x: x[0][0:2], lambda x: matrix(sx0[0:2]).T,[x],x0,name="x[0:2] on dense matrix")
    self.numpyEvaluationCheck(lambda x: x[0][1], lambda x: matrix(sx0[1]).T,[x],x0,name="x[1]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][-1], lambda x: matrix(sx0[-1]).T,[x],x0,name="x[-1]",setx0=[sx0])
    
  def test_getinputExpr(self):
    self.message("outputExpr/inputExpr")
    x=MX("x",2,3)
    f = MXFunction([x],[3*x]) 
    g = MXFunction([f.inputExpr(0)],[6*f.outputExpr(0)]) 
    
    f.init()
    g.init()
    n=[1,2,3,4,5,6]
    f.setInput(n)
    f.evaluate()
    g.setInput(n)
    g.evaluate()
    checkarray(self,6*f.output().toArray(),g.output().toArray(),"slicing(trans)")
    
  def test_scalarMX(self):
      x=MX("x")
      x0=0.738
      self.numpyEvaluationCheckPool(self.pool,[x],x0,name="scalarMX")
      
      self.numpyEvaluationCheckPool(self.matrixpool,[x],x0,name="scalarMX")
      
  def test_MXJacobian(self):
    self.message("MX(1,1) unary operation, jacobian")
    self.Jpool=FunctionPool()
    self.message("SXMatrix(1,1) unary operation, jacobian")
    x=MX("x")
    x0=array([[0.738]])

    def fmod(f,x):
      #f.setOption("ad_mode","forward")
      f.setOption("numeric_jacobian",True)
      f.init()
      J=f.jacobian()
      J.init()
      return J
      
    self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="MX unary operations, jacobian",fmod=fmod)
    
    def fmod(f,x):
      #f.setOption("ad_mode","reverse")
      f.setOption("numeric_jacobian",True)
      f.init()
      J=f.jacobian()
      J.init()
      return J
      
    self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="MX unary operations, jacobian",fmod=fmod)
    
  def test_MXJacobians(self):
      self.message("MX(3,1) unary operation, jacobian")
      x=msym("x",3,1)
      
      x0=array([0.738,0.9,0.3])

      def fmod(f,x):
        #f.setOption("ad_mode","forward")
        f.setOption("numeric_jacobian",True)
        f.init()
        J=f.jacobian()
        J.init()
        return J
        
      self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="MX unary operations, jacobian",fmod=fmod)
      
      def fmod(f,x):
        #f.setOption("ad_mode","reverse")
        f.setOption("numeric_jacobian",True)
        f.init()
        J=f.jacobian()
        J.init()
        return J
        
      self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="MX unary operations, jacobian",fmod=fmod)

  def test_MXbinary(self):
      self.message("MX binary operations")
      x=MX("x",3,2)
      y=MX("x",3,2)
      x0=array([[0.738,0.2],[ 0.1,0.39 ],[0.99,0.999999]])
      y0=array([[1.738,0.6],[ 0.7,12 ],[0,-6]])
      self.numpyEvaluationCheckPool(self.matrixbinarypool,[x,y],[x0,y0],name="MX")

  def test_MXSparse(self):
      self.message("MX unary operations, sparse")
      sp=CRSSparsity(3,4,[1,2,1],[0,2,2,3])
      
      x=MX("x",sp)
      if scipy_available:
        x0=DMatrix(3,4,[1,2,1],[0,2,2,3],[0.738,0.1,0.99]).toCsr_matrix()
        
        self.numpyEvaluationCheckPool(self.pool,[x],array(x0.todense()),name="MX",setx0=x0,excludeflags={'nozero'})
        self.numpyEvaluationCheckPool(self.matrixpool,[x],array(x0.todense()),name="MX",setx0=x0)
      else:
        x0=DMatrix(3,4,[1,2,1],[0,2,2,3],[0.738,0.1,0.99]).toArray()
        
        self.numpyEvaluationCheckPool(self.pool,[x],x0,name="MX",setx0=x0)
        self.numpyEvaluationCheckPool(self.matrixpool,[x],x0,name="MX",setx0=x0)
      
  def test_MXbinarySparse(self):
      self.message("SXMatrix binary operations")
      spx=CRSSparsity(3,4,[1,2,1],[0,2,2,3])
      spy=CRSSparsity(3,4,[0,2,3],[0,2,2,3])
      xx=MX("x",spx)
      yy=MX("y",spy)
      if scipy_available:
        x0=DMatrix(3,4,[1,2,1],[0,2,2,3],[0.738,0.1,0.99]).toCsr_matrix()
        y0=DMatrix(3,4,[0,2,3],[0,2,2,3],[1.738,0.7,-6]).toCsr_matrix()
        
        self.numpyEvaluationCheckPool(self.matrixbinarypool,[xx,yy],[array(x0.todense()),array(y0.todense())],name="MX",setx0=[x0,y0])
      else:
        x0=DMatrix(3,4,[1,2,1],[0,2,2,3],[0.738,0.1,0.99]).toArray()
        y0=DMatrix(3,4,[0,2,3],[0,2,2,3],[1.738,0.7,-6]).toArray()
        
        self.numpyEvaluationCheckPool(self.matrixbinarypool,[xx,yy],[x0,y0],name="MX",setx0=[x0,y0])

  def test_symbolcheck(self):
    self.message("Check if non-symbolic inputs are caught")
    self.assertRaises(RuntimeError, lambda : SXFunction([MX(0)],[MX("x")]))

  def test_unite(self):
    self.message("unite operation")
    import numpy
    numpy.random.seed(42)
    xn = numpy.random.random((3,4))
    x=MX(3,4)
    y=MX("x",3,4)
    z=unite(x,y)
    f = MXFunction([y],[z])
    f.init()
    f.setInput(xn)
    f.evaluate()
    self.checkarray(f.getOutput(),xn,"unite dense")
 
    spx=CRSSparsity(3,4,[1,2,1],[0,2,2,3])
    spy=CRSSparsity(3,4,[0,2,2],[0,1,2,3])

    nx=DMatrix(spx,0)
    for k in range(nx.size()):
      nx[k]= numpy.random.rand()
    ny=DMatrix(spy,0)
    for k in range(nx.size()):
      ny[k]= numpy.random.rand()
      
    nxn = nx.toArray()
    nyn = ny.toArray()
    x=MX("x",spx)
    y=MX("y",spy)
    z=unite(x,y)

    f = MXFunction([x,y],[z])
    f.init()
    f.setInput(nx,0)
    f.setInput(ny,1)
    f.evaluate()
    self.checkarray(f.getOutput(),nxn+nyn,"unite sparse")
     
  def test_imatrix_index(self):
    self.message("IMatrix indexing")
    X = MX("x",2,2)
    Y = X[IMatrix([[0,2],[1,1],[3,3]])]
    
    f = MXFunction([X],[Y])
    f.init()
    f.setInput([1,2,3,4])
    f.evaluate()
    
    self.checkarray(f.getOutput(),array([[1,3],[2,2],[4,4]]),"IMatrix indexing")
    
    Y = X[:,:]
    Y[IMatrix([[0,2]])] = DMatrix([[9,8]])
    
    f = MXFunction([X],[Y])
    f.init()
    f.setInput([1,2,3,4])
    f.evaluate()
    
    self.checkarray(f.getOutput(),array([[9,2],[8,4]]),"IMatrix indexing assignment")
    
  
  def test_IMatrix_index_slice(self):
    self.message("IMatrix combined with slice")

    A = IMatrix(2,2)
    A[0,0] = 0
    A[1,1] = 1
    A[0,1] = 2
    A[1,0] = 0
    
    
    B = MX(DMatrix([[1,2,3],[4,5,6],[7,8,9],[10,11,12]]))
    F = MX(DMatrix([[1,2],[4,5]]))

    f = MXFunction([],[B[:,A]])
    f.init()
    f.evaluate()

    self.checkarray(f.getOutput(),DMatrix([[1,3],[1,2],[4,6],[4,5],[7,9],[7,8],[10,12],[10,11]]),"B[:,A]")
    
    f = MXFunction([],[B[A,:]])
    f.init()
    f.evaluate()
    self.checkarray(f.getOutput(),DMatrix([[1,7,2,8,3,9],[1,4,2,5,3,6]]),"B[A,:]")
    
    self.assertRaises(Exception, lambda : F[:,A])
    
    f = MXFunction([],[B[A,1]])
    f.init()
    f.evaluate()
    
    self.checkarray(f.getOutput(),DMatrix([[2,8],[2,5]]),"B[A,1]")

    f = MXFunction([],[B[1,A]])
    f.init()
    f.evaluate()
    
    self.checkarray(f.getOutput(),DMatrix([[4,6],[4,5]]),"B[1,A]")

    B = MX(DMatrix([[1,2,3],[4,5,6],[7,8,9],[10,11,12]]))
    A = IMatrix([2,0])
    
    B[1,A] = DMatrix([20,21])

    f = MXFunction([],[B])
    f.init()
    f.evaluate()

    self.checkarray(f.getOutput(),DMatrix([[1,2,3],[21,5,20],[7,8,9],[10,11,12]]),"B[1,A] setter")


    B = MX(DMatrix([[1,2,3],[4,5,6],[7,8,9],[10,11,12]]))
    A = IMatrix([2,0])
    
    B[A,1] = DMatrix([20,21])

    f = MXFunction([],[B])
    f.init()
    f.evaluate()

    self.checkarray(f.getOutput(),DMatrix([[1,21,3],[4,5,6],[7,20,9],[10,11,12]]),"B[A,:] setter")
    
    
    B = MX(DMatrix([[1,2,3],[4,5,6],[7,8,9],[10,11,12]]))
    A = IMatrix([2,0])
    
    B[A,:] = DMatrix([[20,21,22],[24,25,26]])

    f = MXFunction([],[B])
    f.init()
    f.evaluate()

    self.checkarray(f.getOutput(),DMatrix([[24,25,26],[4,5,6],[20,21,22],[10,11,12]]),"B[A,:] setter")
    
  def test_IMatrix_IMatrix_index(self):
    self.message("IMatrix IMatrix index")

    A = IMatrix(2,2)
    A[0,0] = 0
    A[1,1] = 1
    A[0,1] = 2
    
    B = IMatrix(2,2)
    B[0,0] = 2
    B[1,1] = 1
    B[0,1] = 0
    
    C = MX(DMatrix([[1,2,3],[4,5,6],[7,8,9],[10,11,12]]))
    F = MX(DMatrix([[1,2],[4,5]]))

    f = MXFunction([],[C[A,B]])
    f.init()
    f.evaluate()
    
    self.checkarray(f.getOutput(),DMatrix([[3,7],[0,5]]),"C[A,B]")
    self.assertRaises(Exception, lambda : F[A,B])

  def test_subsass(self):
     self.message("Check subscripted assignment")
     
     X = MX("x",2,2)
     X[0,0]=MX(5)
     X[0,0]=5
     X[:,0]=8

     x=MX("X",3,4)
     import numpy
     numpy.random.seed(42)
     xn = numpy.random.random((3,4))
     r = numpy.zeros((7,8))
     y=MX.zeros(7,8)
     y[1:4,[2,4,6,7]]=x
     r[1:4,[2,4,6,7]]=xn
     fy = MXFunction([x],[y])
     fy.init()
     fy.setInput(xn)
     fy.evaluate()
     
     self.checkarray(fy.getOutput(),r,"subscripted assigment")
     
     y=MX(7,8)
     y[1:4,[2,4,6,7]]=x
     r[1:4,[2,4,6,7]]=xn
     fy = MXFunction([x],[y])
     fy.init()
     fy.setInput(xn)
     fy.evaluate()
     self.checkarray(fy.getOutput(),r,"subscripted assigment")
     
     
     kl=[2,4,5,8]

     s=y.sparsity()
     for k in kl:
       r[s.getRow()[k],s.col()[k]]=1.0
     
     y[kl]=MX(1)
     fy = MXFunction([x],[y])
     fy.init()
     fy.setInput(xn)
     fy.evaluate()
     self.checkarray(fy.getOutput(),r,"subscripted assigment")
     
     y[kl]=x[[0,1,2,3]]
     s=y.sparsity()
     sx=x.sparsity()
     cnt=0
     for k in kl:
       r[s.getRow()[k],s.col()[k]]=xn[sx.getRow()[cnt],sx.col()[cnt]]
       cnt+=1
     fy = MXFunction([x],[y])
     fy.init()
     fy.setInput(xn)
     fy.evaluate()
     self.checkarray(fy.getOutput(),r,"subscripted assigment")
     
     self.message(":fwdSeed")
     
     x=MX("X",3,1)
     xn = numpy.random.random((3,1))
     y=x**3
     y[1]=x[0]**4
     fy = MXFunction([x],[y])
     fy.init()
     fy.setInput(xn)
     fy.setFwdSeed([1,0,0])
     fy.evaluate(1,0)
     self.checkarray(fy.getOutput(),matrix([xn[0,0]**3,xn[0,0]**4,xn[2,0]**3]).T,"subscripted assigment")
     self.checkarray(fy.getFwdSens(),matrix([3*xn[0,0]**2,4*xn[0,0]**3,0]).T,"subscripted assigment")
     
  def test_erase(self):
    self.message("Erase function")
    self.message(":dense")
    y=MX("Y",7,8)
    import numpy
    r=2*numpy.ones((7,8))
    r[1:4,[2,4,6,7]]=numpy.zeros((3,4))
    z = y *2
    z.erase([1,2,3],[2,4,6,7])
    f = MXFunction([y],[z])
    f.init()
    f.setInput([1]*56)
    e = f.getOutput()
    self.checkarray(f.getOutput(),e,"erase")
    self.message(":sparse")

  def test_MXalgebraDense(self):
    self.message("Test some dense algebraic properties of matrices")
    # issue 96
    n = 3
    m = 4
    import numpy
    numpy.random.seed(42)
    A_ = numpy.random.random((m,n))
    A = MX("A",m,n)
    b_ = numpy.random.random((m,1))
    b = MX("b",m,1)
    C_ = numpy.random.random((m,m))
    C = MX("C",m,m)
    D_ = numpy.random.random((m,n))
    D = MX("D",m,n)
    e_ = numpy.random.random((m,1))
    e = MX("e",m,1)
    x_ = numpy.random.random((n,1))
    x = MX("x",n,1)
    
    Axb = mul(A,x)+b
    Dxe = mul(D,x)+e
    a = mul(mul(trans(Axb),C),Dxe)
    
    f = MXFunction([x,A,b,C,D,e],[a])
    f.init()
    f.setInput(x_,0)
    f.setInput(A_,1)
    f.setInput(b_,2)
    f.setInput(C_,3)
    f.setInput(D_,4)
    f.setInput(e_,5)
    f.evaluate()
    
    f_ = dot(dot((dot(A_,x_)+b_).T,C_),(dot(D_,x_)+e_))
    
    self.checkarray(f.getOutput(),f_,"evaluation")
    
    
    J_ = dot(dot((dot(D_,x_)+e_).T,C_.T),A_) + dot(dot((dot(A_,x_)+b_).T,C_),D_)
    
    for mode in ["forward", "reverse"]:
      #f.setOption("ad_mode",mode)
      f.init()
      J = f.jacobian()
      J.init()
      J.setInput(x_,0)
      J.setInput(A_,1)
      J.setInput(b_,2)
      J.setInput(C_,3)
      J.setInput(D_,4)
      J.setInput(e_,5)
      J.evaluate()
      
      self.checkarray(J.getOutput(),J_,"evaluation")
      
  def test_MXalgebraSparse(self):
    self.message("Test some sparse algebraic properties of matrices")
    if not(scipy_available):
      return
    # issue 97
    n = 3
    m = 4
    import numpy
    numpy.random.seed(42)
    
    def randsparsity(m,n):
      sp = CRSSparsity(m,n)
      for i in range((n*m)/2):
        sp.getNZ(numpy.random.randint(m),numpy.random.randint(n))
      return sp
      
    def gentest(m,n):
      As = randsparsity(m,n)
      A_ = DMatrix(As)
      for k in range(As.size()):
        A_[k]= numpy.random.rand()
      A = MX("A",As)
      return (A_.toCsr_matrix(),A)
    
    (A_,A)=gentest(m,n)
    (b_,b)=gentest(m,1)
    (C_,C)=gentest(m,m)
    (D_,D)=gentest(m,n)
    (e_,e)=gentest(m,1)
    x_ = numpy.random.random((n,1))
    x = MX("x",n,1)
    
    Axb = mul(A,x)+b
    Dxe = mul(D,x)+e
    a = mul(mul(trans(Axb),C),Dxe)
    
    f = MXFunction([x,A,b,C,D,e],[a])
    f.init()
    f.setInput(x_,0)
    f.input(1).set(A_)
    f.input(2).set(b_)
    f.input(3).set(C_)
    f.input(4).set(D_)
    f.input(5).set(e_)
    f.evaluate()


    Axb_ = A_*x_+b_
    Dxe_ = D_*x_+e_
    
    f_ = Axb_.T*C_*Dxe_
    
    self.checkarray(f.getOutput(),f_,"evaluation")
    

    J_ = (D_*x_+e_).T*C_.T*A_ + (A_*x_+b_).T*C_*D_
    
    for mode in ["forward", "reverse"]:
      #f.setOption("ad_mode","forward")
      f.init()
      J = f.jacobian()
      J.init()
      J.setInput(x_,0)
      J.input(1).set(A_)
      J.input(2).set(b_)
      J.input(3).set(C_)
      J.input(4).set(D_)
      J.input(5).set(e_)
      J.evaluate()
      
      self.checkarray(J.getOutput(),J_,"evaluation")

  #@unittest.skipIf(not(scipy_available))
  def test_MXalgebraSparseSparse(self):
    if not(scipy_available):
        return
    self.message("Test some sparse algebraic properties of matrices")
    n = 8
    m = 10
    import numpy
    numpy.random.seed(42)
    
    def randsparsity(m,n):
      sp = CRSSparsity(m,n)
      for k in range((n*m)/2):
        i = numpy.random.randint(m)
        j = numpy.random.randint(n)
        if not(i == m/2):
          if n==1 or not(j == n/2):
            sp.getNZ(i,j)
      return sp
      
    def gentest(m,n):
      As = randsparsity(m,n)
      A_ = DMatrix(As)
      for k in range(As.size()):
        A_[k]= numpy.random.rand()
      A = MX("A",As)
      return (A_.toCsr_matrix(),A)
    
    (A_,A)=gentest(m,n)
    (b_,b)=gentest(m,1)
    (C_,C)=gentest(m,m)
    (D_,D)=gentest(m,n)
    (e_,e)=gentest(m,1)
    x_ = numpy.random.random((n,1))
    x = MX("x",n,1)
    
    Axb = mul(A,x)+b
    Dxe = mul(D,x)+e
    a = mul(mul(trans(Axb),C),Dxe)
    
    f = MXFunction([x,A,b,C,D,e],[a])
    f.init()
    f.setInput(x_,0)
    f.input(1).set(A_)
    f.input(2).set(b_)
    f.input(3).set(C_)
    f.input(4).set(D_)
    f.input(5).set(e_)
    f.evaluate()


    Axb_ = A_*x_+b_
    Dxe_ = D_*x_+e_
    
    f_ = Axb_.T*C_*Dxe_
    
    self.checkarray(f.getOutput(),f_,"evaluation")
    

    J_ = (D_*x_+e_).T*C_.T*A_ + (A_*x_+b_).T*C_*D_
    
    for mode in ["forward", "reverse"]:
      #f.setOption("ad_mode","forward")
      f.init()
      J = f.jacobian()
      J.init()
      J.setInput(x_,0)
      J.input(1).set(A_)
      J.input(2).set(b_)
      J.input(3).set(C_)
      J.input(4).set(D_)
      J.input(5).set(e_)
      J.evaluate()
      
      self.checkarray(J.getOutput(),J_,"evaluation")
      
      
  def test_chaining(self):
    self.message("Chaining SX and MX together")
    x=SX("x")
    y=x**3
    f=SXFunction([x],[y])
    f.init()
    J=f.jacobian()
    J.init()
    
    X=MX("X")
    F=MXFunction([X],J.call([X]))
    F.init()
    
    
    x_=1.7
    F.setInput([x_])
    F.setFwdSeed(1)
    F.setAdjSeed(1)
    F.evaluate(1,1)
    self.checkarray(F.getOutput(),3*x_**2,"Chaining eval")
    self.checkarray(F.getFwdSens(),6*x_,"Chaining fwd")
    self.checkarray(F.getAdjSens(),6*x_,"Chaining adj")
    
  def test_issue104(self):
    self.message("regression test #104")
    x = MX("x")

    F = x**2

    f = MXFunction([x],[F])
    f.init()
    f.setInput([-1])
    f.setFwdSeed([1])
    f.setAdjSeed([1])
    f.evaluate(1,1)
    self.checkarray(f.getFwdSens(),-2,"regression")
    self.checkarray(f.getAdjSens(),-2,"regression")

    f.init()
    f.setInput([0])
    f.setFwdSeed([1])
    f.setAdjSeed([1])
    f.evaluate(1,1)
    self.checkarray(f.getFwdSens(),0,"regression")
    self.checkarray(f.getAdjSens(),0,"regression")
    
    
    x = MX("x",2,1)
    F = x**2
    f = MXFunction([x],[F])
    f.init()
    f.setInput([-1,-1])
    f.setFwdSeed([1,0])
    f.setAdjSeed([1,0])
    f.evaluate(1,1)
    print f.getFwdSens()
    self.checkarray(f.getFwdSens(),matrix([-2,0]).T,"regression")
    self.checkarray(f.getAdjSens(),matrix([-2,0]).T,"regression")

    f.init()
    f.setInput([0,0])
    f.setFwdSeed([1,0])
    f.setAdjSeed([1,0])
    f.evaluate(1,1)
    self.checkarray(f.getFwdSens(),matrix([0,0]).T,"regression")
    self.checkarray(f.getAdjSens(),matrix([0,0]).T,"regression")
    
    x = MX("x")
    y = MX("y")

    F = x**y

    f = MXFunction([x,y],[F])
    f.init()
    f.setInput([-1],0)
    f.setInput([2],1)
    f.setFwdSeed([1])
    f.setAdjSeed([1])
    f.evaluate(1,1)
    self.assertTrue(isnan(f.getFwdSens()[0]))
    self.assertTrue(isnan(f.getAdjSens(1)[0]))
    
    F = constpow(x,y)

    f = MXFunction([x,y],[F])
    f.init()
    f.setInput([-1],0)
    f.setInput([2],1)
    f.setFwdSeed([1])
    f.setAdjSeed([1])
    f.evaluate(1,1)
    self.checkarray(f.getFwdSens(),-2,"regression")
    self.checkarray(f.getAdjSens(),-2,"regression")
    
  def test_issue107(self):
    self.message("Regression test for issue 107: +=")
    x=MX("x")
    y=MX("y")

    z=x
    z+=y
    
    self.assertTrue(isSymbolic(x))
    self.assertFalse(isSymbolic(z))
    
  def test_MXd_trivial(self):
    self.message("symbolic variables and constants jac")
    X =  MX("X",10)
    V =  MX("V")
    f =  MXFunction([X,V],[X,MX.eye(3)])
    f.init()
    self.assertTrue(isinstance(f.jac(0,0),MX))
    self.assertEqual(f.jac(0,0).size(),10)
    self.assertEqual(f.jac(0,0).size1(),10)
    self.assertEqual(f.jac(0,0).size2(),10)
    
    g = MXFunction([],[f.jac(0,0)]);g.init();g.evaluate()
    self.checkarray(g.getOutput(),eye(10),"unit matrix")
    
    g = MXFunction([],[f.jac(0,1)]);g.init();g.evaluate()
    self.checkarray(g.getOutput(),zeros((9,10)),"zero matrix")
    
    g = MXFunction([],[f.jac(1,0)]);g.init();g.evaluate()
    self.checkarray(g.getOutput(),zeros((10,1)),"zero matrix")
    
    g = MXFunction([],[f.jac(1,1)]);g.init();g.evaluate()
    self.checkarray(g.getOutput(),zeros((9,1)),"zero matrix")
    
  def test_MXd_substractionl(self):
    self.message("substraction jac")
    V =  MX("V")
    X =  MX("X")
    f =  MXFunction([X,V],[X-V])
    f.init()
    
    g = MXFunction([],[f.jac(0,0)]);g.init();g.evaluate()
    self.checkarray(g.getOutput(),ones((1,1)),"one")

    g = MXFunction([],[f.jac(1,0)]);g.init();g.evaluate()
    self.checkarray(g.getOutput(),-ones((1,1)),"one")
    
    f =  MXFunction([X,V],[V-X])
    f.init()
    
    g = MXFunction([],[f.jac(0,0)]);g.init();g.evaluate()
    self.checkarray(g.getOutput(),-ones((1,1)),"one")

    g = MXFunction([],[f.jac(1,0)]);g.init();g.evaluate()
    self.checkarray(g.getOutput(),ones((1,1)),"one")
    
  def test_MXd_mapping(self):
    self.message("mapping jac")
    X = MX("X",3)
    Y = MX("Y",2)
    f = MXFunction([X,Y],[vertcat([X,Y])])
    f.init()
    J = f.jac(0,0)
    JJ = DMatrix(J.sparsity(),1)
    self.checkarray(JJ,vstack((eye(3),zeros((2,3)))),"diag")
    J = f.jac(1,0)
    JJ = DMatrix(J.sparsity(),1)
    self.checkarray(JJ,vstack((zeros((3,2)),eye(2))),"diag")
    
    
  def test_MatrixAlgebraTableDense(self):
    self.message("Table of derivatives https://ccrma.stanford.edu/~dattorro/matrixcalc.pdf")
    
    n = m = K = L = k = 3
    t_ = 0.3
    mu_ = 0.13

    def gentest(m,n):
      A = MX("A",m,n)
      return (DMatrix(numpy.random.random((m,n))),A)
 
    (a_,a) = gentest(m,1)
    (b_,b) = gentest(m,1)
    (x_,x) = gentest(n,1)
    (y_,y) = gentest(n,1)
    (A_,A) = gentest(m,n)
    (B_,B) = gentest(m,n)
    (C_,C) = gentest(m,n)
    (X_,X) = gentest(K,L)
    (Y_,Y) = gentest(K,L)
    t = MX(t_)
    mu = MX(mu_)
    
    ins = [a,b,x,y,A,B,C,X,Y]
    ins_ = [a_,b_,x_,y_,A_,B_,C_,X_,Y_]
    
    def grad(y,x):
      f = MXFunction(ins,[x])
      f.init()
      J = f.jacobian([i for i in range(len(ins)) if ins[i] is y][0])
      J.init()
      if x.shape[0]==1 and x.shape[1]==1:
        return (J.call(ins)[0].T).reshape(y.shape)
      return J.call(ins)[0].T

    def eye(n):
      return DMatrix(numpy.eye(n))
    
    Axb = mul(A,x)-b
    ab = mul(a,b.T)
    tests = [
    (grad(x,x),eye(k)),
    (grad(x,x.T),eye(k)),
    (grad(x,Axb),A.T),
    #(grad(x,Axb.T),A)   incorrect?
    (grad(x,mul(Axb.T,Axb)),2*mul(A.T,Axb)),
    #(grad(x,norm_2(Axb)),mul(A.T,Axb)/norm_2(Axb)), #  norm_2 not implemented
    (grad(x,mul(mul(x.T,A),x)+2*mul(mul(x.T,B),y)+mul(mul(y.T,C),y)),mul((A+A.T),x)+2*mul(B,y)),
    #(grad(x,mul(a.T,mul(x.T,x)*b)),2*mul(mul(x,a.T),b))
    (grad(X,X),eye(k**2)),
    #(grad(X,X.T),eye(k**2))
    (grad(X,mul(a.T,mul(X,b))),ab),
    (grad(X,mul(b.T,mul(X.T,a))),ab),
    (grad(X,mul(a.T,mul(mul(X,X),b))),mul(X.T,ab)+mul(ab,X.T)),
    (grad(X,mul(a.T,mul(mul(X.T,X),b))),mul(X,ab + ab.T)),
    (grad(x,x*mu),MX(eye(k))*mu),
    (grad(X,c.trace(X*mu)),MX(eye(k))*mu),
    (grad(X,c.trace(mul(X.T,Y))),Y),
    (grad(X,c.trace(mul(Y,X.T))),Y),
    (grad(X,c.trace(mul(Y.T,X))),Y),
    (grad(X,c.trace(mul(X,Y.T))),Y),
    (grad(X,c.trace(mul(a.T,mul(X,b)))),ab)
    #(grad(X,log(c.det(X))),c.inv(X_)),
    ]

    cnt = 0
    for symbol, solution in tests:
      f = MXFunction(ins,[symbol])
      f.init()
      for i in range(len(ins_)):
        f.setInput(ins_[i],i)
      f.evaluate()
      g = MXFunction(ins,[solution])
      g.init()
      for i in range(len(ins_)):
        g.setInput(ins_[i],i)
      g.evaluate()
      self.checkarray(f.getOutput(),g.getOutput(),"#%d" % cnt )
      cnt+=1
    
  def test_issue134(self):
    self.message("Test issue #134")

    x = MX("x",2,2)
    y = MX(2,2)
    
    x_ = DMatrix([[1,2],[3,4]])

    f = MXFunction([x],[x+y])
    f.init()
    f.setInput(x_,0)
    f.evaluate(0,0) # this should not throw a segfault
    self.checkarray(f.getOutput(),x_,"issue 134")
    f.evaluate(1,1) # this should not throw a segfault

    f = MXFunction([x],[y+x])
    f.init()
    f.setInput(x_,0)
    f.evaluate(0,0) # this should not throw a segfault
    self.checkarray(f.getOutput(),x_,"issue 134")
    f.evaluate(1,1) # this should not throw a segfault
    
    x = MX("x",1,1)
    y = MX(1,1)

    x_ = 7.1
    
    f = MXFunction([x],[x+y])
    f.init()
    f.setInput(x_,0)
    f.evaluate(0,0) # this should not throw a segfault
    self.checkarray(f.getOutput(),x_,"issue 134")
    f.evaluate(1,1) # this should not throw a segfault
    
    f = MXFunction([x],[y+x])
    f.init()
    f.setInput(x_,0)
    f.evaluate(0,0) # this should not throw a segfault
    self.checkarray(f.getOutput(),x_,"issue 134")
    f.evaluate(1,1) # this should not throw a segfault
        
  # 2-norms currently not supported
  #def test_Norm2(self):
    #self.message("Norm_2")
    #X=MX("x",5,1)
    
    #nums = matrix([1,2,3,0,-1]).T

    #F =MXFunction([X],[norm_2(X)])

    #J = F.jacobian(0,0)
    #J.init()

    #J.setInput(nums)
    #J.evaluate()
    #self.checkarray(J.getOutput(),nums.T/linalg.norm(nums),"Norm_2")

    #J = Jacobian(F,0,0)
    #J.setOption("ad_mode","reverse")
    #J.init()

    #J.setInput(nums)
    #J.evaluate()
    
    #self.checkarray(J.getOutput(),nums.T/linalg.norm(nums),"Norm_2")


    #J = MXFunction([X],[F.jac(0)[0]])
    #J.init()

    #J.setInput(nums)
    #J.evaluate()
    
    #self.checkarray(J.getOutput(),nums.T/linalg.norm(nums),"Norm_2")
    
  # Removed since a normed squared is not a norm
  #def test_Norm22(self):
    #self.message("Norm_22")
    #X=MX("x",5,1)
    
    #nums = matrix([1,2,3,0,-1]).T

    #F =MXFunction([X],[norm_22(X)])

    #J = Jacobian(F,0,0)
    #J.setOption("ad_mode","forward")
    #J.init()

    #J.setInput(nums)
    #J.evaluate()
    #self.checkarray(J.getOutput(),2*nums.T,"Norm_22 fwd")

    #J = Jacobian(F,0,0)
    #J.setOption("ad_mode","reverse")
    #J.init()

    #J.setInput(nums)
    #J.evaluate()
    
    #self.checkarray(J.getOutput(),2*nums.T,"Norm_22 adj")
        
    #J = MXFunction([X],[F.jac(0)[0]])
    #J.init()

    #J.setInput(nums)
    #J.evaluate()
    
    #self.checkarray(J.getOutput(),2*nums.T,"Norm_22 jac")

  # 1-norms currently not supported
  #def test_Norm1(self):
    #self.message("Norm_1")
    #X=MX("x",3,1)
    
    #nums = matrix([6,-3,0]).T

    #F =MXFunction([X],[norm_1(X)])

    #J = Jacobian(F,0,0)
    #J.setOption("ad_mode","forward")
    #J.init()

    #J.setInput(nums)
    #J.evaluate()
    #self.checkarray(J.getOutput(),matrix([1,-1,nan]),"Norm_1")

    #J = Jacobian(F,0,0)
    #J.setOption("ad_mode","reverse")
    #J.init()

    #J.setInput(nums)
    #J.evaluate()
    
    #self.checkarray(J.getOutput(),matrix([1,-1,nan]),"Norm_1")
    
  def test_null(self):
    self.message("MXFunction null")
    x = msym("x")

    f = MXFunction([x],[x**2,MX()])
    f.init()

    self.assertEqual(f.output(1).shape[0],0)
    self.assertEqual(f.output(1).shape[1],0)
    f.evaluate(1,1)
    self.assertEqual(f.fwdSens(1).shape[0],0)
    self.assertEqual(f.fwdSens(1).shape[1],0)
    self.assertEqual(f.adjSeed(1).shape[0],0)
    self.assertEqual(f.adjSeed(1).shape[1],0)
    
    f = MXFunction([x,MX()],[x**2,MX()])
    f.init()

    self.assertEqual(f.output(1).shape[0],0)
    self.assertEqual(f.output(1).shape[1],0)
    f.evaluate(1,1)
    self.assertEqual(f.fwdSens(1).shape[0],0)
    self.assertEqual(f.fwdSens(1).shape[1],0)
    self.assertEqual(f.adjSeed(1).shape[0],0)
    self.assertEqual(f.adjSeed(1).shape[1],0)
    
    r = f.call([x,MX()])
    self.assertTrue(r[1].isNull())

    r = f.call([MX(),MX()])
    self.assertTrue(r[1].isNull())
    
    #self.assertRaises(Exception,lambda : f.eval([x,x]))
    #self.assertRaises(Exception,lambda : f.eval([[],[]]))
    
  def test_issue184(self):
    self.message("Regression test issue #184")
    x = MX("x", 3)
    y = x[0:0]
    self.assertEqual(y.size(),0)

  def test_indexingOutOfBounds(self):
    self.message("Indexing out of bounds")
    y = DMatrix.zeros(4, 5) 
    self.assertRaises(RuntimeError,lambda : y[12,0] )
    self.assertRaises(RuntimeError,lambda : y[12,12] )
    self.assertRaises(RuntimeError,lambda : y[0,12] )
    self.assertRaises(RuntimeError,lambda : y[12,:] )
    self.assertRaises(RuntimeError,lambda : y[12:15,0] )
    self.assertRaises(RuntimeError,lambda : y[:,12] )
    self.assertRaises(RuntimeError,lambda : y[0,12:15] )
    y[-1,2]
    self.assertRaises(RuntimeError,lambda : y[-12,2] )
    y[-3:-1,2]
    self.assertRaises(RuntimeError,lambda : y[-12:-9,2] )
    
    def test():
      y[12,0] = 0
    self.assertRaises(RuntimeError,test)
    def test():
      y[12,12] = 0
    self.assertRaises(RuntimeError,test)
    def test():
      y[0,12] = 0
    self.assertRaises(RuntimeError,test)
    def test():
      y[12,:] = 0
    self.assertRaises(RuntimeError,test)
    def test():
      y[12:15,0] = 0
    self.assertRaises(RuntimeError,test)
    def test():
      y[:,12] = 0
    self.assertRaises(RuntimeError,test)
    def test():
      y[0,12:15] = 0
    self.assertRaises(RuntimeError,test)
    y[-1,2] = 0
    def test():
      y[-12,2] = 0
    self.assertRaises(RuntimeError,test)
    y[-3:-1,2] = 0
    def test():
      y[-12:-9,2] = 0
    self.assertRaises(RuntimeError,test)
    
  def test_indexinglimits(self):
    self.message("Limits of indexing")
    y = casadi.MX("y", 3) 
    self.assertRaises(RuntimeError,lambda : y[[0, 5]] )
    try:
      y[[0, 5]] = MX("a")
      self.assertTrue(False)
    except RuntimeError:
      pass
    y[[0, 2]]
    y[[0, 2]] = MX("a")
    
  def test_printLimiting(self):
    self.message("printLimiting")

    x = MX("x")
    for i in range(100):
      x = sin(x)*x
      
    self.assertTrue(len(str(x)) <  4*MX.getMaxNumCallsInPrint())

    MX.setMaxNumCallsInPrint(5)
    self.assertTrue(len(str(x)) <  100)

    MX.setMaxNumCallsInPrint()

  def test_mul(self):
    A = MX(DMatrix.ones((4,3)))
    B = MX(DMatrix.ones((3,8)))
    C = MX(DMatrix.ones((8,7)))
    
    self.assertRaises(RuntimeError,lambda : mul([]))
    
    D = mul([A])
    
    self.assertEqual(D.shape[0],4)
    self.assertEqual(D.shape[1],3)

    D = mul([A,B])
    
    self.assertEqual(D.shape[0],4)
    self.assertEqual(D.shape[1],8)
    
    D = mul([A,B,C])
    
    self.assertEqual(D.shape[0],4)
    self.assertEqual(D.shape[1],7)
 
  def test_truth(self):
    self.message("Truth values")
    self.assertRaises(Exception, lambda : bool(msym("x")))
    #self.assertRaises(Exception, lambda : bool(msym("x")>0))
    self.assertTrue(bool(MX(1)))
    self.assertFalse(bool(MX(0)))
    self.assertTrue(bool(MX(0.2)))
    self.assertTrue(bool(MX(-0.2)))
    self.assertRaises(Exception, lambda : bool(MX(DMatrix([2.0,3]))))
    self.assertRaises(Exception, lambda : bool(MX()))


  def test_MXbool(self):
    self.message("bool")
    
    xy = MX("x",2)
    x = xy[0]
    y = xy[1]
    
    f = MXFunction([xy],[vertcat([logic_and(x,y),logic_or(x,y),logic_not(x)])])
    f.init()
    
    
    for t1 in [0,1]:
      for t2 in [0,1]:
        T1 = t1!=0
        T2 = t2!=0
        f.setInput([t1,t2])
        f.evaluate()
        self.checkarray(f.getOutput(),DMatrix([T1 and T2,T1 or T2,not T1]),"bool(%d,%d): %s" % (t1,t2,str(f.getOutput())))

  def test_MXineq(self):
    self.message("SX ineq")
    
    xy = MX("x",2)
    x = xy[0]
    y = xy[1]
    
    
    f = MXFunction([xy],[vertcat([x<y,x<=y,x>=y,x==y,x!=y])])
    f.init()
    
    
    for t1 in [-10,0.1,0,1,10]:
      for t2 in [-10,0.1,0,1,10]:
        T1 = t1
        T2 = t2
        f.setInput([t1,t2])
        f.evaluate()
        self.checkarray(f.getOutput(),DMatrix([T1 < T2,T1 <= T2, T1 >= T2, T1 == T2, T1 != T2]),"ineq(%d,%d)" % (t1,t2))

  def test_if_else_zero(self):
    x = MX("x")
    y = if_else(x,5,0)
    f = MXFunction([x],[y])
    f.init()
    f.setInput(1)
    f.evaluate()
    self.assertTrue(f.getOutput()==5,"if_else_zero %s " % str(f.getOutput()))
    f.setInput(0)
    f.evaluate()
    self.assertTrue(f.getOutput()==0,"if_else_zero")
    
    
  def test_if_else(self):
    x = MX("x")
    y = if_else(x,1,2)
    f = MXFunction([x],[y])
    f.init()
    f.setInput(1)
    f.evaluate()
    self.assertTrue(f.getOutput()==1,"if_else")
    f.setInput(0)
    f.evaluate()
    self.assertTrue(f.getOutput()==2,"if_else")
    
    # Check sensitivities
    
    x0 = 2.1
    dx = 0.3
    y = if_else(x>1,x**2,x**3)
    f = MXFunction([x],[y])
    f.init()
    f.setInput(x0)
    f.setFwdSeed(dx)
    f.setAdjSeed(dx)
    f.evaluate(1,1)
    self.checkarray(f.getOutput(),x0**2,"if_else sens")
    self.checkarray(f.getFwdSens(),2*x0*dx,"if_else sens")
    self.checkarray(f.getAdjSens(),2*x0*dx,"if_else sens")
    
    x0 = -2.1
    dx = 0.3
    
    f.setInput(x0)
    f.setFwdSeed(dx)
    f.setAdjSeed(dx)
    f.evaluate(1,1)
    self.checkarray(f.getOutput(),x0**3,"if_else sens")
    self.checkarray(f.getFwdSens(),3*(-x0)**2*dx,"if_else sens")
    self.checkarray(f.getAdjSens(),3*(-x0)**2*dx,"if_else sens")
    
  def test_regression491(self):
    self.message("regression #491")
    u = ssym("u")
    x = ssym("x")

    F = SXFunction([u,x],[u+1/x])
    F.init()

    U = MX("U")

    X = F.call([U,U])[0]
    G = F.call([U,X])[0]

    for kk in range(3):
      gfcn = 0
      if kk==0:
        tmp = MXFunction([U],[G])
        tmp.init()
        gfcn = tmp.expand()
      else:
        gfcn = MXFunction([U],[G])
        gfcn.setOption("numeric_jacobian",kk==1)
      gfcn.setOption("ad_mode","reverse")
      gfcn.init()
      J = gfcn.jacobian()
      J.init()
      J.setInput(1)
      J.evaluate()
      self.assertAlmostEqual(J.getOutput(),1,9)

  def test_ticket(self):
    J = [] + msym("x")
    J = msym("x") + []
        
  def test_jacobian_tools(self):
    self.message("jacobian")
    
    X = msym("X")

    Y = jacobian(X**2,X)
    
    f = MXFunction([X],[Y])
    f.init()
    
    f.setInput(2.3)
    f.evaluate()
    
    self.assertAlmostEqual(f.getOutput(),4.6)
    
  def test_reshape(self):
    self.message("reshape")
    X = msym("X",10)

    i = IMatrix(sp_tril(3),range(6))

    i.printDense()
    print vecNZ(i)

    T = X[i]

    f = MXFunction([X],[vecNZ(T)**2])
    f.init()
    f.setInput(range(10))
    f.evaluate()
    
    self.checkarray(IMatrix([0,1,9,4,16,25]),f.getOutput())

    Y = msym("Y",10)

    ff = MXFunction([Y],f.evalMX([Y]))
    ff.init()
    ff.setInput(range(10))
    ff.evaluate()

    self.checkarray(IMatrix([0,1,9,4,16,25]),ff.getOutput())
    
    J = MXFunction([X],[f.jac()])
    J.init()
    J.setInput(range(10))
    J.evaluate()
    
    i = horzcat([diag([0,2,4,6,8,10]),IMatrix.zeros(6,4)])
    i[[2,3],:] = i[[3,2],:]
    
    self.checkarray(i,J.getOutput())
    
    f = MXFunction([X],[vecNZ(T)**2])
    f.setOption("ad_mode","reverse")
    f.init()
    
    J = MXFunction([X],[f.jac()])
    J.init()
    J.setInput(range(10))
    J.evaluate()
    
    i = horzcat([diag([0,2,4,6,8,10]),IMatrix.zeros(6,4)])
    i[[2,3],:] = i[[3,2],:]
    
    self.checkarray(i,J.getOutput())
    
  def test_vertcat(self):
    self.message("vertcat")
    X = msym("X",10)

    T = vertcat([X[4],X[2]])

    f = MXFunction([X],[T**2])
    f.init()
    f.setInput(range(10))
    f.evaluate()
    
    self.checkarray(IMatrix([16,4]),f.getOutput())

    Y = msym("Y",10)

    ff = MXFunction([Y],f.evalMX([Y]))
    ff.init()
    ff.setInput(range(10))
    ff.evaluate()

    self.checkarray(IMatrix([16,4]),ff.getOutput())
    
    J = MXFunction([X],[f.jac()])
    J.init()
    J.setInput(range(10))
    J.evaluate()
    
    i = IMatrix.zeros(2,10)
    i[0,4] = 8
    i[1,2] = 4
    
    self.checkarray(i,J.getOutput())
    
    f = MXFunction([X],[T**2])
    f.setOption("ad_mode","reverse")
    f.init()
    
    J = MXFunction([X],[f.jac()])
    J.init()
    J.setInput(range(10))
    J.evaluate()
    
    self.checkarray(i,J.getOutput())
  
  def test_blockcat(self):
    x = msym("x")
    
    y = blockcat([[x,2*x],[3*x,4*x]])
    f = MXFunction([x],[y])
    f.init()
    f.setInput(3)
    f.evaluate()
    self.checkarray(f.getOutput(),DMatrix([[3,6],[9,12]]))
    
    
  def test_veccats(self):
    x= msym("x",2)
    self.assertTrue(hash(vec(x))==hash(x))
    self.assertTrue(hash(vecNZ(x))==hash(x))
    
  def test_constmxmul(self):
    0.1*MX.ones(2)

  def test_isRegular(self):
    self.assertTrue(isRegular(MX(DMatrix([0,1]))))
    self.assertFalse(isRegular(MX(DMatrix([0,Inf]))))
    with self.assertRaises(Exception):
      self.assertFalse(isRegular(msym("x",2)))

  def test_blkdiag(self):
    C = blkdiag([MX(DMatrix(([[-1.4,-3.2],[-3.2,-28]]))),DMatrix([[15,-12,2.1],[-12,16,-3.8],[2.1,-3.8,15]]),1.8,-4.0])
    self.assertTrue(isinstance(C,MX))
    r = DMatrix([[-1.4,-3.2,0,0,0,0,0],[-3.2,-28,0,0,0,0,0],[0,0,15,-12,2.1,0,0],[0,0,-12,16,-3.8,0,0],[0,0,2.1,-3.8,15,0,0],[0,0,0,0,0,1.8,0],[0,0,0,0,0,0,-4]])
    makeSparse(r)
    f = MXFunction([],[C])
    f.init()
    f.evaluate()
    
    self.checkarray(f.getOutput(),r)
    
  def test_tril2symm(self):
    x = msym("x",sp_tril(3))
    f = MXFunction([x],[tril2symm(x)])
    f.init()
    f.setInput(range(6))
    f.evaluate()
    self.checkarray(f.getOutput(),DMatrix([[0,1,3],[1,2,4],[3,4,5]]))
    
  def test_sparsity_indexing(self):
    self.message("sparsity")

    B_ = DMatrix([[1,2,3,4,5],[6,7,8,9,10]])
    B = msym("B",2,5)
    
    A = IMatrix([[1,1,0,0,0],[0,0,1,0,0]])
    makeSparse(A)
    sp = A.sparsity()
    import copy
    
    def meval(m):
      f = MXFunction([B],[m])
      f.init()
      f.setInput(B_)
      f.evaluate()
      return f.getOutput() 
    
    self.checkarray(meval(B[sp]),DMatrix([[1,2,0,0,0],[0,0,8,0,0]]),"sparsity indexing")

    Bmod = copy.copy(B)
    Bmod[sp] = -4
    
    self.checkarray(meval(Bmod),DMatrix([[-4,-4,3,4,5],[6,7,-4,9,10]]),"sparsity indexing assignement")

    Bmod = copy.copy(B)
    Bmod[sp] = 2*B
    
    self.checkarray(meval(Bmod),DMatrix([[2,4,3,4,5],[6,7,16,9,10]]),"Imatrix indexing assignement")
    
    self.assertRaises(Exception, lambda : B[sp_dense(4,4)])

  def test_getSymbols(self):
    a = msym("a")
    b = msym("b")
    c = msym("c")
    e = cos(a*b) + c
    w = getSymbols(e)
    self.assertEqual(len(w),3)
    if CasadiOptions.getSimplificationOnTheFly():
      self.assertTrue(isEqual(w[0],a))
      self.assertTrue(isEqual(w[1],b))
      self.assertTrue(isEqual(w[2],c))
    
  def test_iter(self):
    self.assertEqual(len(list(msym("x",2))),2)

  def test_vertcat_empty(self):
    a = MX(DMatrix(0,2))
    v = vertcat([a,a])
    
    self.assertEqual(v.size1(),0)
    self.assertEqual(v.size2(),2)
    
    a = MX(DMatrix(2,0))
    v = vertcat([a,a])
    
    self.assertEqual(v.size1(),4)
    self.assertEqual(v.size2(),0)
    
if __name__ == '__main__':
    unittest.main()
