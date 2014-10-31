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
    x = MX.sym("x",1,3)
    z=vertcat([x*(i+1) for i in range(8)])
    f = MXFunction([x],[ztf(z)])
    f.init()
    L=[1,2,3]
    f.setInput(L,0)
    f.evaluate()
    zt = f.getOutput(0).toArray()
    zr = array([[L[0]*(i+1),L[1]*(i+1),L[2]*(i+1)] for i in range(8)])
    checkarray(self,zrf(zr),zt,name)
    return (zt,zrf(zr))

def checkMXoperations2(self,ztf,zrf,name):
    x = MX.sym("x",3,1)
    z = horzcat([x*i for i in range(8)])
    f = MXFunction([x],[ztf(z)])
    f.init()
    L=[1,2,3]
    f.setInput(L,0)
    f.evaluate()
    zt = f.getOutput(0).toArray()
    zr = array([[L[0]*i,L[1]*i,L[2]*i] for i in range(8)]).T
    checkarray(self,zrf(zr),zt,name)
    return zt

def checkMXoperations3(self,ztf,zrf,name):
    x = MX.sym("x",3,1)
    p = horzcat([x[0,0],x[1,0],x[2,0]])
    z = vertcat([p*i for i in range(8)])
    f = MXFunction([x],[ztf(z)])
    f.init()
    L=[1,2,3]
    f.setInput(L,0)
    f.evaluate()
    zt = f.getOutput(0).toArray()
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
    self.pool.append(lambda x: log(x[0]),log,"log",flags={'nozero'})
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
    self.matrixbinarypool.append(lambda a: fmax(a[0],a[1]),lambda a: fmax(a[0],a[1]),"fmin")

    self.matrixbinarypool.append(lambda a: fmin(a[0],a[1]),lambda a: fmin(a[0],a[1]),"fmax")
    self.matrixbinarypool.append(lambda a: mul(a[0],a[1].T),lambda a: dot(a[0],a[1].T),"mul(Matrix,Matrix.T)")
    self.matrixbinarypool.append(lambda a: arctan2(a[0],a[1]),lambda a: arctan2(a[0],a[1]),"arctan2")
    #self.matrixbinarypool.append(lambda a: inner_mul(a[0],trans(a[1])),lambda a: dot(a[0].T,a[1]),name="inner_mul(Matrix,Matrix)") 
    self.matrixbinarypool.append(lambda a: mul(a[0],a[1].T),lambda a: dot(a[0],a[1].T),"mul(Matrix,Matrix.T)")
    
  def test_MX1(self):
    self.message("MX constructor")
    x = MX.sym("x",2,3)
    self.assertEqual(x.size1(),2,"MX fails to indicate its size1")
    self.assertEqual(x.size2(),3,"MX fails to indicate its size2")

  def test_MXvertcat(self):
    self.message("MX vertcat")
    x = MX.sym("x",1,3)
    y = MX.sym("y",1,3)
    z=vertcat((x,y))
    self.assertEqual(z.size1(),2,"MX fails to indicate its size1")
    self.assertEqual(z.size2(),3,"MX fails to indicate its size2")

  def test_MXFunction1(self):
    self.message("MXFunction single input, single output")
    # check if x->2*x
    # evaluates correctly for x=3
    x = MX.sym("x")
    y = 2*x
    f = MXFunction([x],[y])
    self.assertEqual(f.getNumInputs(),1,"MXFunction fails to indicate correct number of inputs")
    self.assertEqual(f.getNumOutputs(),1,"MXFunction fails to indicate correct number of outputs")
    f.init()
    f.setInput(3,0);
    f.evaluate()
    yt = tuple(f.getOutput().data())
    self.assertEqual(type(yt),TupleType,"Output of MXFunction is expected to be tuple of floats")
    self.assertEqual(len(yt),1,"Output of MXFunction was tuple of floats, as expected, but length is incorrect.")
    y=yt[0]
    self.assertEqual(type(y),float,"Output of MXFunction is expected to be tuple of floats")
    self.assertAlmostEqual(y, 2*3,10)

  def test_MXfunction2(self):
    self.message("MXFunction multi input, multi output")
      # check if [x,y]->[y+x,y*x]
    # evaluates correctly for x=3,y=7
    x = MX.sym("x")
    y = MX.sym("y")
    f = MXFunction([x,y],[x+y,y*x])
    self.assertEqual(f.getNumInputs(),2,"MXFunction fails to indicate correct number of inputs")
    self.assertEqual(f.getNumOutputs(),2,"MXFunction fails to indicate correct number of outputs")

    f.init()
    f.setInput(3,0);
    f.setInput(7,1);
    f.evaluate()
    zt1 = tuple(f.getOutput(0).data())
    zt2 = tuple(f.getOutput(1).data())
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
    xy = MX.sym("xy",2)
    f = MXFunction([xy],[xy[0]+xy[1],xy[0]*xy[1]])
    self.assertEqual(f.getNumInputs(),1,"MXFunction fails to indicate correct number of inputs")
    self.assertEqual(f.getNumOutputs(),2,"MXFunction fails to indicate correct number of outputs")
    f.init()
    f.setInput([3,7],0);
    f.evaluate()
    zt1 = tuple(f.getOutput(0).data())
    zt2 = tuple(f.getOutput(1).data())
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
    xy = MX.sym("xy",1,2)
    f = MXFunction([xy],[xy[0,0]+xy[0,1],xy[0,0]*xy[0,1]])
    self.assertEqual(f.getNumInputs(),1,"MXFunction fails to indicate correct number of inputs")
    self.assertEqual(f.getNumOutputs(),2,"MXFunction fails to indicate correct number of outputs")
    f.init()
    f.setInput([3,7],0);
    f.evaluate()
    zt1 = f.getOutput(0).toArray()
    zt2 = f.getOutput(1).toArray()
    
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
    xy = MX.sym("xy",2)
    z=vertcat([xy[0]+xy[1],xy[0]*xy[1]])
    f = MXFunction([xy],[z])
    self.assertEqual(f.getNumInputs(),1,"MXFunction fails to indicate correct number of inputs")
    self.assertEqual(f.getNumOutputs(),1,"MXFunction fails to indicate correct number of outputs")
    f.init()
    f.setInput([3,7],0);
    f.evaluate()
    zt=f.getOutput(0).toArray()
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
    xy = MX.sym("xy",2)
    z=horzcat([xy[0]+xy[1],xy[0]*xy[1]])
    f = MXFunction([xy],[z])
    self.assertEqual(f.getNumInputs(),1,"MXFunction fails to indicate correct number of inputs")
    self.assertEqual(f.getNumOutputs(),1,"MXFunction fails to indicate correct number of outputs")
    f.init()
    f.setInput([3,7],0);
    f.evaluate()
    zt = f.getOutput(0).toArray()
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
    x=MX.sym("x")
    y=MX.sym("y")

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
    x = SXElement.sym("x")
    f = SXFunction([x],[x])
    f.init()
    f.setInput([3],0)
    f.evaluate()
    self.assertAlmostEqual(f.getOutput(0)[0,0], 3,10)

  def test_identityMX(self):
    self.message("identity MXFunction")
    x = MX.sym("x")
    f  = MXFunction([x],[x])
    f.init()
    f.setInput([3],0)
    f.evaluate()
    self.assertAlmostEqual(f.getOutput(0)[0,0], 3,10)
    
  def test_MXorder(self):
    self.message("MXFunction order of non-zero elements")
    x = MX.sym("x",2,3)
    f = MXFunction([x],[x+x])

    self.assertEqual(f.getNumInputs(),1,"MXFunction fails to indicate correct number of inputs")
    self.assertEqual(f.getNumOutputs(),1,"MXFunction fails to indicate correct number of outputs")
    f.init()
    L=[1,2,3,4,5,6]
    f.setInput(L,0)
    f.evaluate()
    zt = f.getOutput(0).toArray()
    self.assertEqual(zt.shape[0],2,"Output of MXFunction is of wrong shape.")
    self.assertEqual(zt.shape[1],3,"Output of MXFunction is of wrong shape.")
      
    Lr=reshape(L,(2,3),'F')
    for i in range(2):
      for j in range(3):
        self.assertAlmostEqual(Lr[i,j]*2, zt[i,j],10)
    
  def test_trans(self):
    self.message("trans")
    a = MX.sparse(0,1)
    b = a.T
    self.assertEquals(b.size1(),1)
    self.assertEquals(b.size2(),0)
    
  def test_MXtrans(self):
    self.message("trans(MX)")
    x = MX.sym("x",2,3)
    z=x.T
    self.assertEqual(z.size1(),3,"Vec returns MX of wrong dimension")
    self.assertEqual(z.size2(),2,"Vec returns MX of wrong dimension")
    f = MXFunction([x],[z])
    self.assertEqual(f.getNumInputs(),1,"MXFunction fails to indicate correct number of inputs")
    self.assertEqual(f.getNumOutputs(),1,"MXFunction fails to indicate correct number of outputs")
    f.init()
    L=[1,2,3,4,5,6]
    f.setInput(L,0)
    f.evaluate()
    zt = f.getOutput(0).toArray()
    
    ztr=reshape(zt,(3,2))
    Lr=reshape(L,(2,3),'F')
    for i in range(2):
      for j in range(3):
        self.assertAlmostEqual(Lr[i,j], ztr[j,i],10)
      
  def test_MXvec(self):

    u = DMatrix([[10*j+i for i in range(3)] for j in range(4) ])

    U = MX.sym("u",u.shape)

    f = MXFunction([U],[vec(U)])
    f.init()
    f.setInput(u)
    f.evaluate()
    
    self.checkarray(vec(u),f.getOutput(),"vec")
    
  def test_MXvecNZ(self):

    u = DMatrix.sparse(4,3)
    u[0,0] = 7
    u[1,1] = 8
    u[2,2] = 6
    u[1,0] = 9
    u[0,1] = 11
    
    U = MX.sym("u",u.sparsity())

    f = MXFunction([U],[vecNZ(U)])
    f.init()
    f.setInput(u)
    f.evaluate()
    
    self.checkarray(vecNZ(u),f.getOutput(),"vec")

  def test_MXreshape(self):
    self.message("reshape(MX)")
    x = MX.sym("x",2,3)
    z=c.reshape(x,(1,6))
    self.assertEqual(z.size1(),1,"Vec returns MX of wrong dimension")
    self.assertEqual(z.size2(),6,"Vec returns MX of wrong dimension")
    f = MXFunction([x],[z])
    self.assertEqual(f.getNumInputs(),1,"MXFunction fails to indicate correct number of inputs")
    self.assertEqual(f.getNumOutputs(),1,"MXFunction fails to indicate correct number of outputs")
    f.init()
    L=[1,2,3,4,5,6]
    f.setInput(L,0)
    f.evaluate()
    zt = f.getOutput(0).toArray()
    for i in range(len(L)):
      self.assertAlmostEqual(L[i], zt[0,i],10)
  
  def test_MXcompose(self):
    self.message("compositions of vec, trans, reshape with vertcat")
    checkMXoperations(self,lambda x: x,lambda x: x,'vertcat')
    checkMXoperations(self,lambda x: x.T,lambda x: x.T,'trans(vertcat)')
    checkMXoperations(self,lambda x: x.T.T,lambda x: x,'trans(trans(vertcat))')
    checkMXoperations(self,lambda x: vec(x.T),lambda x: reshape(x,(prod(x.shape),1)),'vec(trans(vertcat))')
    checkMXoperations(self,lambda x: vec(x).T,lambda x: reshape(x.T,(prod(x.shape),1)).T,'vec(trans(vertcat))')
    checkMXoperations(self,lambda x: c.reshape(x.T,(6,4)).T,lambda x: reshape(x,(4,6)),'reshape(vertcat)')
    checkMXoperations(self,lambda x: c.reshape(x,(6,4)).T,lambda x: reshape(x.T,(4,6)),'reshape(trans(vertcat))') 
    checkMXoperations(self,lambda x: c.reshape(x.T,(6,4)),lambda x: reshape(x,(4,6)).T,'trans(reshape(vertcat))') 

  def test_MXcompose2(self):
    self.message("compositions of vec, trans, reshape with horzcat")
    checkMXoperations2(self,lambda x: x,lambda x: x,'horzcat')
    checkMXoperations2(self,lambda x: x.T,lambda x: x.T,'trans(horzcat)')
    checkMXoperations2(self,lambda x: x.T.T,lambda x: x,'trans(trans(horzcat))')
    checkMXoperations2(self,lambda x: vec(x.T),lambda x: reshape(x,(prod(x.shape),1)),'vec(trans(horzcat))')
    checkMXoperations2(self,lambda x: vec(x).T,lambda x: reshape(x.T,(prod(x.shape),1)).T,'vec(trans(horzcat))')
    checkMXoperations2(self,lambda x: c.reshape(x.T,(6,4)).T,lambda x: reshape(x,(4,6)),'reshape(horzcat)')
    checkMXoperations2(self,lambda x: c.reshape(x,(6,4)).T,lambda x: reshape(x.T,(4,6)),'reshape(trans(horzcat))') 
    checkMXoperations2(self,lambda x: c.reshape(x.T,(6,4)),lambda x: reshape(x,(4,6)).T,'trans(reshape(horzcat))') 

  def test_MXcompose3(self):
    self.message("compositions of vec, trans, reshape with vertcat")
    checkMXoperations3(self,lambda x: x,lambda x: x,'snippet')
    checkMXoperations3(self,lambda x: x.T,lambda x: x.T,'trans(snippet)')
    checkMXoperations3(self,lambda x: x.T.T,lambda x: x,'trans(trans(snippet))')
    checkMXoperations3(self,lambda x: vec(x.T),lambda x: reshape(x,(prod(x.shape),1)),'vec(trans(snippet))')
    checkMXoperations3(self,lambda x: vec(x).T,lambda x: reshape(x.T,(prod(x.shape),1)).T,'vec(trans(snippet))')
    checkMXoperations3(self,lambda x: c.reshape(x.T,(6,4)).T,lambda x: reshape(x,(4,6)),'reshape(snippet)')
    checkMXoperations3(self,lambda x: c.reshape(x,(6,4)).T,lambda x: reshape(x.T,(4,6)),'reshape(trans(snippet))') 
    checkMXoperations3(self,lambda x: c.reshape(x.T,(6,4)),lambda x: reshape(x,(4,6)).T,'trans(reshape(snippet))') 

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
    x = MX.sym("x",3,2)
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
    self.numpyEvaluationCheck(lambda x: x[0].nz[[0,2,3]], lambda x: matrix([x[0,0],x[2,0],x[0,1]]).T,[x],x0,name="x[[0,2,3]]")

    self.numpyEvaluationCheck(lambda x: x[0].nz[1], lambda x: matrix(x.T.ravel()[1]).T,[x],x0,name="x[1] on dense matrix")
    self.numpyEvaluationCheck(lambda x: x[0].nz[-1], lambda x: matrix(x.ravel()[-1]).T,[x],x0,name="x[-1] on dense matrix")
    self.numpyEvaluationCheck(lambda x: x[0][[0,1],0:1],lambda x: x[[0,1],0:1],[x],x0,name='x[:,0:1]')
    self.numpyEvaluationCheck(lambda x: x[0][0:1,[0,1]],lambda x: x[0:1,[0,1]],[x],x0,name='x[0:1,:]')
    
    self.message(":sparse")
      
    sp=Sparsity(4,3,[0,2,2,3],[1,2,1])
    x=MX.sym("X",sp)
    sx0=[0.738,0.39,0.99]
    x0=DMatrix(Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.39,0.99]).toArray()
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
    self.numpyEvaluationCheck(lambda x: x[0].nz[[2,1]], lambda x: matrix([x[1,2],x[2,0]]).T,[x],x0,name="x[[2,1]]")
    self.numpyEvaluationCheck(lambda x: x[0].nz[0:2], lambda x: matrix(sx0[0:2]).T,[x],x0,name="x[0:2] on dense matrix")
    self.numpyEvaluationCheck(lambda x: x[0].nz[1], lambda x: matrix(sx0[1]).T,[x],x0,name="x[1]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0].nz[-1], lambda x: matrix(sx0[-1]).T,[x],x0,name="x[-1]",setx0=[sx0])
    
  def test_getinputExpr(self):
    self.message("outputExpr/inputExpr")
    x=MX.sym("x",2,3)
    f = MXFunction([x],[3*x]) 
    g = MXFunction([f.inputExpr(0)],[6*f.outputExpr(0)]) 
    
    f.init()
    g.init()
    n=[1,2,3,4,5,6]
    f.setInput(n)
    f.evaluate()
    g.setInput(n)
    g.evaluate()
    checkarray(self,6*f.getOutput().toArray(),g.getOutput().toArray(),"slicing(trans)")
    
  def test_scalarMX(self):
      x=MX.sym("x")
      x0=0.738
      self.numpyEvaluationCheckPool(self.pool,[x],x0,name="scalarMX")
      
      self.numpyEvaluationCheckPool(self.matrixpool,[x],x0,name="scalarMX")
      
  def test_MXJacobian(self):
    self.message("MX(1,1) unary operation, jacobian")
    self.Jpool=FunctionPool()
    self.message("SX(1,1) unary operation, jacobian")
    x=MX.sym("x")
    x0=array([[0.738]])

    def fmod(f,x):
      #f.setOption("ad_mode","forward")
      f.init()
      J=f.jacobian()
      J.init()
      return J
      
    self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="MX unary operations, jacobian",fmod=fmod)
    
    def fmod(f,x):
      #f.setOption("ad_mode","reverse")
      f.init()
      J=f.jacobian()
      J.init()
      return J
      
    self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="MX unary operations, jacobian",fmod=fmod)
    
  def test_MXJacobians(self):
      self.message("MX(3,1) unary operation, jacobian")
      x=MX.sym("x",3,1)
      
      x0=array([0.738,0.9,0.3])

      def fmod(f,x):
        #f.setOption("ad_mode","forward")
        f.init()
        J=f.jacobian()
        J.init()
        return J
        
      self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="MX unary operations, jacobian",fmod=fmod)
      
      def fmod(f,x):
        #f.setOption("ad_mode","reverse")
        f.init()
        J=f.jacobian()
        J.init()
        return J
        
      self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="MX unary operations, jacobian",fmod=fmod)

  def test_MXbinary(self):
      self.message("MX binary operations")
      x=MX.sym("x",3,2)
      y=MX.sym("x",3,2)
      x0=array([[0.738,0.2],[ 0.1,0.39 ],[0.99,0.999999]])
      y0=array([[1.738,0.6],[ 0.7,12 ],[0,-6]])
      self.numpyEvaluationCheckPool(self.matrixbinarypool,[x,y],[x0,y0],name="MX")

  def test_MXSparse(self):
      self.message("MX unary operations, sparse")
      sp=Sparsity(4,3,[0,2,2,3],[1,2,1])
      
      x=MX.sym("x",sp)
      if scipy_available:
        x0=DMatrix(Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).toCsc_matrix()
        
        self.numpyEvaluationCheckPool(self.pool,[x],array(x0.todense()),name="MX",setx0=x0,excludeflags={'nozero'})
        self.numpyEvaluationCheckPool(self.matrixpool,[x],array(x0.todense()),name="MX",setx0=x0)
      else:
        x0=DMatrix(Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).toArray()
        
        self.numpyEvaluationCheckPool(self.pool,[x],x0,name="MX",setx0=x0,excludeflags={'nozero'})
        self.numpyEvaluationCheckPool(self.matrixpool,[x],x0,name="MX",setx0=x0)
      
  def test_MXbinarySparse(self):
      self.message("SX binary operations")
      spx=Sparsity(4,3,[0,2,2,3],[1,2,1])
      spy=Sparsity(4,3,[0,2,2,3],[0,2,3])
      xx=MX.sym("x",spx)
      yy=MX.sym("y",spy)
      if scipy_available:
        x0=DMatrix(Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).toCsc_matrix()
        y0=DMatrix(Sparsity(4,3,[0,2,2,3],[0,2,3]),[1.738,0.7,-6]).toCsc_matrix()
        
        self.numpyEvaluationCheckPool(self.matrixbinarypool,[xx,yy],[array(x0.todense()),array(y0.todense())],name="MX",setx0=[x0,y0])
      else:
        x0=DMatrix(Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).toArray()
        y0=DMatrix(Sparsity(4,3,[0,2,2,3],[0,2,3]),[1.738,0.7,-6]).toArray()
        
        self.numpyEvaluationCheckPool(self.matrixbinarypool,[xx,yy],[x0,y0],name="MX",setx0=[x0,y0])

  def test_symbolcheck(self):
    self.message("Check if non-symbolic inputs are caught")
    self.assertRaises(RuntimeError, lambda : SXFunction([MX(0)],[MX.sym("x")]))

  def test_unite(self):
    self.message("unite operation")
    import numpy
    numpy.random.seed(42)
    xn = numpy.random.random((3,4))
    x=MX.sparse(3,4)
    y=MX.sym("x",3,4)
    z=unite(x,y)
    f = MXFunction([y],[z])
    f.init()
    f.setInput(xn)
    f.evaluate()
    self.checkarray(f.getOutput(),xn,"unite dense")
 
    spx=Sparsity(4,3,[0,2,2,3],[1,2,1])
    spy=Sparsity(4,3,[0,1,2,3],[0,2,2])

    nx=DMatrix(spx,0)
    for k in range(nx.size()):
      nx.nz[k]= numpy.random.rand()
    ny=DMatrix(spy,0)
    for k in range(nx.size()):
      ny.nz[k]= numpy.random.rand()
      
    nxn = nx.toArray()
    nyn = ny.toArray()
    x=MX.sym("x",spx)
    y=MX.sym("y",spy)
    z=unite(x,y)

    f = MXFunction([x,y],[z])
    f.init()
    f.setInput(nx,0)
    f.setInput(ny,1)
    f.evaluate()
    self.checkarray(f.getOutput(),nxn+nyn,"unite sparse")
     
  def test_imatrix_index(self):
    self.message("IMatrix indexing")
    X = MX.sym("x",2,2)
    Y = X.nz[IMatrix([[0,2],[1,1],[3,3]])]
    
    f = MXFunction([X],[Y])
    f.init()
    f.setInput([1,2,3,4])
    f.evaluate()
    
    self.checkarray(f.getOutput(),array([[1,3],[2,2],[4,4]]),"IMatrix indexing")
    
    Y = X[:,:]
    Y.nz[IMatrix([[0,2]])] = DMatrix([[9,8]])
    
    f = MXFunction([X],[Y])
    f.init()
    f.setInput([1,2,3,4])
    f.evaluate()
    
    self.checkarray(f.getOutput(),array([[9,8],[2,4]]),"IMatrix indexing assignment")
    
  
  def test_IMatrix_index_slice(self):
    self.message("IMatrix combined with slice")

    A = IMatrix.sparse(2,2)
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

    A = IMatrix.sparse(2,2)
    A[0,0] = 0
    A[1,1] = 1
    A[0,1] = 2
    
    B = IMatrix.sparse(2,2)
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
     
     X = MX.sym("x",2,2)
     X[0,0]=MX(5)
     X[0,0]=5
     X[:,0]=8

     x=MX.sym("X",3,4)
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
     
     y=MX.sparse(7,8)
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
       r[s.row()[k],s.getCol()[k]]=1.0
     
     y.nz[kl]=MX(1)
     fy = MXFunction([x],[y])
     fy.init()
     fy.setInput(xn)
     fy.evaluate()
     self.checkarray(fy.getOutput(),r,"subscripted assigment")
     
     y.nz[kl]=x.nz[[0,1,2,3]]
     s=y.sparsity()
     sx=x.sparsity()
     cnt=0
     for k in kl:
       r[s.row()[k],s.getCol()[k]]=xn[sx.row()[cnt],sx.getCol()[cnt]]
       cnt+=1
     fy = MXFunction([x],[y])
     fy.init()
     fy.setInput(xn)
     fy.evaluate()
     self.checkarray(fy.getOutput(),r,"subscripted assigment")
          
  def test_erase(self):
    self.message("Erase function")
    self.message(":dense")
    y=MX.sym("Y",7,8)
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
    A = MX.sym("A",m,n)
    b_ = numpy.random.random((m,1))
    b = MX.sym("b",m,1)
    C_ = numpy.random.random((m,m))
    C = MX.sym("C",m,m)
    D_ = numpy.random.random((m,n))
    D = MX.sym("D",m,n)
    e_ = numpy.random.random((m,1))
    e = MX.sym("e",m,1)
    x_ = numpy.random.random((n,1))
    x = MX.sym("x",n,1)
    
    Axb = mul(A,x)+b
    Dxe = mul(D,x)+e
    a = mul(mul(Axb.T,C),Dxe)
    
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
      sp = Sparsity.sparse(m,n)
      for i in range((n*m)/2):
        sp.getNZ(numpy.random.randint(m),numpy.random.randint(n))
      return sp
      
    def gentest(m,n):
      As = randsparsity(m,n)
      A_ = DMatrix(As)
      for k in range(As.size()):
        A_.nz[k]= numpy.random.rand()
      A = MX.sym("A",As)
      return (A_.toCsc_matrix(),A)
    
    (A_,A)=gentest(m,n)
    (b_,b)=gentest(m,1)
    (C_,C)=gentest(m,m)
    (D_,D)=gentest(m,n)
    (e_,e)=gentest(m,1)
    x_ = numpy.random.random((n,1))
    x = MX.sym("x",n,1)
    
    Axb = mul(A,x)+b
    Dxe = mul(D,x)+e
    a = mul(mul(Axb.T,C),Dxe)
    
    f = MXFunction([x,A,b,C,D,e],[a])
    f.init()
    f.setInput(x_,0)
    f.setInput(A_,1)
    f.setInput(b_,2)
    f.setInput(C_,3)
    f.setInput(D_,4)
    f.setInput(e_,5)
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
      J.setInput(A_,1)
      J.setInput(b_,2)
      J.setInput(C_,3)
      J.setInput(D_,4)
      J.setInput(e_,5)
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
      sp = Sparsity.sparse(m,n)
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
        A_.nz[k]= numpy.random.rand()
      A = MX.sym("A",As)
      return (A_.toCsc_matrix(),A)
    
    (A_,A)=gentest(m,n)
    (b_,b)=gentest(m,1)
    (C_,C)=gentest(m,m)
    (D_,D)=gentest(m,n)
    (e_,e)=gentest(m,1)
    x_ = numpy.random.random((n,1))
    x = MX.sym("x",n,1)
    
    Axb = mul(A,x)+b
    Dxe = mul(D,x)+e
    a = mul(mul(Axb.T,C),Dxe)
    
    f = MXFunction([x,A,b,C,D,e],[a])
    f.init()
    f.setInput(x_,0)
    f.setInput(A_,1)
    f.setInput(b_,2)
    f.setInput(C_,3)
    f.setInput(D_,4)
    f.setInput(e_,5)
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
      J.setInput(A_,1)
      J.setInput(b_,2)
      J.setInput(C_,3)
      J.setInput(D_,4)
      J.setInput(e_,5)
      J.evaluate()
      
      self.checkarray(J.getOutput(),J_,"evaluation")
      
      
  def test_chaining(self):
    self.message("Chaining SXElement and MX together")
    x=SXElement.sym("x")
    y=x**3
    f=SXFunction([x],[y])
    f.init()
    J=f.jacobian()
    J.init()
    
    X=MX.sym("X")
    F=MXFunction([X],J.call([X]))
    F.init()
    
    
    x_=1.7
    F.setInput([x_])
    F.evaluate()
    self.checkarray(F.getOutput(),3*x_**2,"Chaining eval")
    
  def test_issue107(self):
    self.message("Regression test for issue 107: +=")
    x=MX.sym("x")
    y=MX.sym("y")

    z=x
    z+=y
    
    self.assertTrue(x.isSymbolic())
    self.assertFalse(z.isSymbolic())
    
  def test_MXd_trivial(self):
    self.message("symbolic variables and constants jac")
    X =  MX.sym("X",10)
    V =  MX.sym("V")
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
    V =  MX.sym("V")
    X =  MX.sym("X")
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
    X = MX.sym("X",3)
    Y = MX.sym("Y",2)
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
      A = MX.sym("A",m,n)
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
            
  # 2-norms currently not supported
  #def test_Norm2(self):
    #self.message("Norm_2")
    #X=MX.sym("x",5,1)
    
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
    
  # 1-norms currently not supported
  #def test_Norm1(self):
    #self.message("Norm_1")
    #X=MX.sym("x",3,1)
    
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
    x = MX.sym("x")

    f = MXFunction([x],[x**2,MX()])
    f.init()

    self.assertEqual(f.getOutput(1).shape[0],0)
    self.assertEqual(f.getOutput(1).shape[1],0)
    f.evaluate()
    
    f = MXFunction([x,MX()],[x**2,MX()])
    f.init()

    self.assertEqual(f.getOutput(1).shape[0],0)
    self.assertEqual(f.getOutput(1).shape[1],0)
    f.evaluate()
    
    r = f.call([x,MX()])
    self.assertTrue(r[1].isEmpty(True))

    r = f.call([MX(),MX()])
    self.assertTrue(r[1].isEmpty(True))
    
    #self.assertRaises(Exception,lambda : f.call([x,x],True))
    #self.assertRaises(Exception,lambda : f.call([[],[]],True))
    
  def test_issue184(self):
    self.message("Regression test issue #184")
    x = MX.sym("x", 3)
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
    y = casadi.MX.sym("y", 3) 
    self.assertRaises(RuntimeError,lambda : y[[0, 5]] )
    try:
      y[[0, 5]] = MX.sym("a")
      self.assertTrue(False)
    except RuntimeError:
      pass
    y[[0, 2]]
    y[[0, 2]] = MX.sym("a")
    
  def test_printLimiting(self):
    self.message("printLimiting")

    x = MX.sym("x")
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
    self.assertRaises(Exception, lambda : bool(MX.sym("x")))
    #self.assertRaises(Exception, lambda : bool(MX.sym("x")>0))
    self.assertTrue(bool(MX(1)))
    self.assertFalse(bool(MX(0)))
    self.assertTrue(bool(MX(0.2)))
    self.assertTrue(bool(MX(-0.2)))
    self.assertRaises(Exception, lambda : bool(MX(DMatrix([2.0,3]))))
    self.assertRaises(Exception, lambda : bool(MX()))


  def test_MXbool(self):
    self.message("bool")
    
    xy = MX.sym("x",2)
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
    
    xy = MX.sym("x",2)
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
    x = MX.sym("x")
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
    x = MX.sym("x")
    y = if_else(x,1,2)
    f = MXFunction([x],[y])
    f.init()
    f.setInput(1)
    f.evaluate()
    self.assertTrue(f.getOutput()==1,"if_else")
    f.setInput(0)
    f.evaluate()
    self.assertTrue(f.getOutput()==2,"if_else")
        
  def test_regression491(self):
    self.message("regression #491")
    u = SX.sym("u")
    x = SX.sym("x")

    F = SXFunction([u,x],[u+1/x])
    F.init()

    U = MX.sym("U")

    X = F.call([U,U])[0]
    G = F.call([U,X])[0]

    for kk in range(2):
      gfcn = 0
      if kk==0:
        tmp = MXFunction([U],[G])
        tmp.init()
        gfcn = tmp.expand()
      else:
        gfcn = MXFunction([U],[G])
      gfcn.setOption("ad_mode","reverse")
      gfcn.init()
      J = gfcn.jacobian()
      J.init()
      J.setInput(1)
      J.evaluate()
      self.assertAlmostEqual(J.getOutput(),1,9)

  def test_ticket(self):
    J = [] + MX.sym("x")
    J = MX.sym("x") + []
        
  def test_jacobian_tools(self):
    self.message("jacobian")
    
    X = MX.sym("X")

    Y = jacobian(X**2,X)
    
    f = MXFunction([X],[Y])
    f.init()
    
    f.setInput(2.3)
    f.evaluate()
    
    self.assertAlmostEqual(f.getOutput(),4.6)
    
  def test_reshape(self):
    self.message("reshape")
    X = MX.sym("X",10)

    i = IMatrix(Sparsity.tril(3),range(6))

    i.printDense()
    print vecNZ(i.T)

    T = X[i]

    f = MXFunction([X],[vecNZ(T.T)**2])
    f.init()
    f.setInput(range(10))
    f.evaluate()
    
    self.checkarray(IMatrix([0,1,9,4,16,25]),f.getOutput())

    Y = MX.sym("Y",10)

    ff = MXFunction([Y],f.call([Y],True))
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
    
    f = MXFunction([X],[vecNZ(T.T)**2])
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
    X = MX.sym("X",10)

    T = vertcat([X[4],X[2]])

    f = MXFunction([X],[T**2])
    f.init()
    f.setInput(range(10))
    f.evaluate()
    
    self.checkarray(IMatrix([16,4]),f.getOutput())

    Y = MX.sym("Y",10)

    ff = MXFunction([Y],f.call([Y],True))
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
    x = MX.sym("x")
    
    y = blockcat([[x,2*x],[3*x,4*x]])
    f = MXFunction([x],[y])
    f.init()
    f.setInput(3)
    f.evaluate()
    self.checkarray(f.getOutput(),DMatrix([[3,6],[9,12]]))
    
    
  def test_veccats(self):
    x= MX.sym("x",2)
    self.assertTrue(hash(vec(x))==hash(x))
    self.assertTrue(hash(vecNZ(x))==hash(x))
    
  def test_constmxmul(self):
    0.1*MX.ones(2)

  def test_isRegular(self):
    self.assertTrue(MX(DMatrix([0,1])).isRegular())
    self.assertFalse(MX(DMatrix([0,Inf])).isRegular())
    with self.assertRaises(Exception):
      self.assertFalse(MX.sym("x",2).isRegular())

  def test_blkdiag(self):
    C = blkdiag([MX(DMatrix(([[-1.4,-3.2],[-3.2,-28]]))),DMatrix([[15,-12,2.1],[-12,16,-3.8],[2.1,-3.8,15]]),1.8,-4.0])
    self.assertTrue(isinstance(C,MX))
    r = DMatrix([[-1.4,-3.2,0,0,0,0,0],[-3.2,-28,0,0,0,0,0],[0,0,15,-12,2.1,0,0],[0,0,-12,16,-3.8,0,0],[0,0,2.1,-3.8,15,0,0],[0,0,0,0,0,1.8,0],[0,0,0,0,0,0,-4]])
    r = sparse(r)
    f = MXFunction([],[C])
    f.init()
    f.evaluate()
    
    self.checkarray(f.getOutput(),r)
    
  def test_tril2symm(self):
    x = MX.sym("x",Sparsity.tril(3))
    f = MXFunction([x],[tril2symm(x)])
    f.init()
    f.setInput(range(6))
    f.evaluate()
    self.checkarray(f.getOutput(),DMatrix([[0,1,2],[1,3,4],[2,4,5]]))
    
  def test_sparsity_indexing(self):
    self.message("sparsity")

    B_ = DMatrix([[1,2,3,4,5],[6,7,8,9,10]])
    B = MX.sym("B",2,5)
    
    A = IMatrix([[1,1,0,0,0],[0,0,1,0,0]])
    A = sparse(A)
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
    
    self.assertRaises(Exception, lambda : B[Sparsity.dense(4,4)])

  def test_getSymbols(self):
    a = MX.sym("a")
    b = MX.sym("b")
    c = MX.sym("c")
    e = cos(a*b) + c
    w = getSymbols(e)
    self.assertEqual(len(w),3)
    if CasadiOptions.getSimplificationOnTheFly():
      self.assertTrue(isEqual(w[0],a))
      self.assertTrue(isEqual(w[1],b))
      self.assertTrue(isEqual(w[2],c))
    
  def test_iter(self):
    self.assertEqual(len(list(MX.sym("x",2))),2)

  @known_bug()
  def test_vertcat_empty(self):
    a = MX(DMatrix(0,2))
    v = vertcat([a,a])
    
    self.assertEqual(v.size1(),0)
    self.assertEqual(v.size2(),2)
    
    a = MX(DMatrix(2,0))
    v = vertcat([a,a])
    
    self.assertEqual(v.size1(),4)
    self.assertEqual(v.size2(),0)
    
  def test_jacobian_empty(self):
    x = MX.sym("x",3)

    s = jacobian(DMatrix.sparse(0,0),x).shape
    self.assertEqual(s[0],0)
    self.assertEqual(s[1],3)

    s = jacobian(x,MX.sym("x",0,4)).shape
    self.assertEqual(s[0],3)
    self.assertEqual(s[1],0)
    
  def test_mul_sparsity(self):

    N = 10
    x = MX.sym("x",N,N)
    y = MX.sym("y",N,N)

    x_ = self.randDMatrix(N,N)
    y_ = self.randDMatrix(N,N)

    filt = Sparsity.diag(N)+Sparsity.triplet(N,N,[1],[3])

    f = MXFunction([x,y],[mul(x,y)])
    f.init()
    f.setInput(x_,0)
    f.setInput(y_,1)
    g = MXFunction([x,y],[mul(x,y,filt)])
    g.init()
    g.setInput(x_,0)
    g.setInput(y_,1)
    
    f.evaluate()
    g.evaluate()
    
    self.checkarray(IMatrix(filt,1),IMatrix(g.getOutput().sparsity(),1))
    
    self.checkarray(f.getOutput()[filt],g.getOutput())
    
  def test_mul_zero_wrong(self):
    with self.assertRaises(RuntimeError):
      mul(MX.sym("X",4,5),MX.zeros(3,2))
      
  def test_vertsplit(self):
    a = MX.sym("X",Sparsity.tril(5))
    v = vertsplit(a,[0,2,4,5])
    
    f = MXFunction([a],v)
    f.init()
    f.setInput(range(5*6/2))

    f.evaluate()
    v = [f.getOutput(i) for i in range(len(v))]
    
    self.assertEqual(len(v),3)
    self.checkarray(v[0],DMatrix([[0,0,0,0,0],[1,5,0,0,0]]))
    self.checkarray(v[1],DMatrix([[2,6,9,0,0],[3,7,10,12,0]]))
    self.checkarray(v[2],DMatrix([[4,8,11,13,14]]))
    
    v = vertsplit(a)
    
    f = MXFunction([a],v)
    f.init()
    f.setInput(range(5*6/2))

    f.evaluate()
    v = [f.getOutput(i) for i in range(len(v))]
    
    self.assertEqual(len(v),a.size1())
    self.checkarray(v[0],DMatrix([[0,0,0,0,0]]))
    self.checkarray(v[1],DMatrix([[1,5,0,0,0]]))
    self.checkarray(v[2],DMatrix([[2,6,9,0,0]]))
    self.checkarray(v[3],DMatrix([[3,7,10,12,0]]))
    self.checkarray(v[4],DMatrix([[4,8,11,13,14]]))
    
    v = vertsplit(a,2)
    
    f = MXFunction([a],v)
    f.init()
    f.setInput(range(5*6/2))

    f.evaluate()
    v = [f.getOutput(i) for i in range(len(v))]
    
    self.assertEqual(len(v),3)
    self.checkarray(v[0],DMatrix([[0,0,0,0,0],[1,5,0,0,0]]))
    self.checkarray(v[1],DMatrix([[2,6,9,0,0],[3,7,10,12,0]]))
    self.checkarray(v[2],DMatrix([[4,8,11,13,14]]))
    
    v = vertsplit(a,[0,0,3,a.size1()])
    
    f = MXFunction([a],v)
    f.init()
    f.setInput(range(5*6/2))

    f.evaluate()
    V = [f.getOutput(i) for i in range(len(v))]
    
    self.assertEqual(len(v),3)
    self.assertEqual(v[0].size1(),0)
    self.assertEqual(v[0].size2(),5)  # why not 5?
    self.checkarray(V[1],DMatrix([[0,0,0,0,0],[1,5,0,0,0],[2,6,9,0,0]]))
    self.checkarray(V[2],DMatrix([[3,7,10,12,0],[4,8,11,13,14]]))
 
  def test_horzsplit(self):
    a = MX.sym("X",Sparsity.tril(5))
    v = horzsplit(a,[0,2,4,5])
    
    f = MXFunction([a],v)
    f.init()
    f.setInput(range(5*6/2))

    f.evaluate()
    v = [f.getOutput(i) for i in range(len(v))]
    self.assertEqual(len(v),3)
    self.checkarray(v[0],DMatrix([[0,0],[1,5],[2,6],[3,7],[4,8]]))
    self.checkarray(v[1],DMatrix([[0,0],[0,0],[9,0],[10,12],[11,13]]))
    self.checkarray(v[2],DMatrix([[0],[0],[0],[0],[14]]))
    
    v = horzsplit(a)
    
    f = MXFunction([a],v)
    f.init()
    f.setInput(range(5*6/2))

    f.evaluate()
    v = [f.getOutput(i) for i in range(len(v))]
    self.assertEqual(len(v),a.size1())
    self.checkarray(v[0],DMatrix([0,1,2,3,4]))
    self.checkarray(v[1],DMatrix([0,5,6,7,8]))
    self.checkarray(v[2],DMatrix([0,0,9,10,11]))
    self.checkarray(v[3],DMatrix([0,0,0,12,13]))
    self.checkarray(v[4],DMatrix([0,0,0,0,14]))
    
    v = horzsplit(a,2)

    f = MXFunction([a],v)
    f.init()
    f.setInput(range(5*6/2))

    f.evaluate()
    v = [f.getOutput(i) for i in range(len(v))]
    
    self.assertEqual(len(v),3)
    self.checkarray(v[0],DMatrix([[0,0],[1,5],[2,6],[3,7],[4,8]]))
    self.checkarray(v[1],DMatrix([[0,0],[0,0],[9,0],[10,12],[11,13]]))
    self.checkarray(v[2],DMatrix([[0],[0],[0],[0],[14]]))
    
    v = horzsplit(a,[0,0,3,a.size2()])
    f = MXFunction([a],v)
    f.init()
    f.setInput(range(5*6/2))

    f.evaluate()
    V = [f.getOutput(i) for i in range(len(v))]
    
    self.assertEqual(len(v),3)
    self.assertEqual(v[0].size1(),5)
    self.assertEqual(v[0].size2(),0)
    self.checkarray(V[1],DMatrix([[0,0,0],[1,5,0],[2,6,9],[3,7,10],[4,8,11]]))
    self.checkarray(V[2],DMatrix([[0,0],[0,0],[0,0],[12,0],[13,14]]))
    
  def test_blocksplit(self):
    a = MX.sym("X",Sparsity.tril(5))
    v = blocksplit(a,[0,2,4,5],[0,1,3,5])
    
    fs = [MXFunction([a],vr) for vr in v]
    for f in fs:
      f.init()
      f.setInput(range(5*6/2))

      f.evaluate()
    v = [[fs[i].getOutput(j) for j in range(3)] for i in range(3)]
    
    self.checkarray(v[0][0],DMatrix([0,1]))
    self.checkarray(v[0][1],DMatrix([[0,0],[5,0]]))
    self.checkarray(v[1][0],DMatrix([2,3]))
    self.checkarray(blockcat(v),f.getInput())

  def test_mxnulloutput(self):
     a = MX.sparse(5,0)
     b = MX.sym("x",2)
     
     f = MXFunction([b],[a])
     f.init()
     c = f.call([b])[0]

     self.assertEqual(c.size1(),5)
     self.assertEqual(c.size2(),0)

     c = f.call([b],True)[0]

     self.assertEqual(c.size1(),5)
     self.assertEqual(c.size2(),0)
     
     a = MX.sparse(0,0)
     b = MX.sym("x",2)
     
     f = MXFunction([b],[a])
     f.init()
     c = f.call([b])[0]

     self.assertEqual(c.size1(),0)
     self.assertEqual(c.size2(),0)

     c = f.call([b],True)[0]

     self.assertEqual(c.size1(),0)
     self.assertEqual(c.size2(),0)
     
  def test_mxnull(self):
     a = MX.sparse(5,0)
     b = MX.sparse(0,3)
     
     c = mul(a,b)
     
     self.assertEqual(c.size(),0)
     
     a = MX.sparse(5,3)
     b = MX.sparse(3,4)
     
     c = mul(a,b)
     
     self.assertEqual(c.size(),0)
     
  def  test_mxnullop(self):
    c = MX.sparse(0,0)
    x = MX.sym("x",2,3)
    
    with self.assertRaises(RuntimeError):
      d = x + c
 
    with self.assertRaises(RuntimeError):
      d = x / c
      
  @slow()
  @memory_heavy()
  def test_MX_shapes(self):
      self.message("MX unary operations")
      
      #self.checkarray(DMatrix(Sparsity.tril(4),1),DMatrix(Sparsity.dense(4,4),1))
      
      for sp in [Sparsity.dense(0,0),Sparsity.dense(0,2),Sparsity.dense(2,0),Sparsity.dense(1,1),Sparsity.dense(2,2), Sparsity(4,3,[0,2,2,3],[1,2,1])]:
        for v in [0,1,0.2]:
          x_ = DMatrix(sp,v)
          
          xx = MX.sym("x",sp.size1(),sp.size2())
          x=xx[sp]
          
          for (casadiop, numpyop,name, flags) in self.pool.zip():
            if 'nozero' in flags and v==0: continue
            r = casadiop([x])
            f = MXFunction([xx],[r])
            f.init()
            f.setInput(v,0)
            f.evaluate()
            
            self.checkarray(f.getOutput(),numpyop(x_))
            
            a = IMatrix(f.getOutput().sparsity(),1)
            b = IMatrix(DMatrix(numpyop(x_)).sparsity(),1)
            
            c = b-a
            if c.size()>0:
              # At least as sparse as DMatrix calculus
              self.assertTrue(min(c)>=0,str([sp,v,name]))

      for sp in [Sparsity.sparse(1,1),Sparsity.dense(1,1),Sparsity.sparse(3,4),Sparsity.dense(3,4), Sparsity(4,3,[0,2,2,3],[1,2,1]).T]:
        for v1 in [0,1,0.2,-0.2]:
          x1_ = DMatrix(sp,v1)
          xx1 = MX.sym("x",sp.size1(),sp.size2())
          x1=xx1[sp]
          xx1s = SX.sym("x",sp.size1(),sp.size2())
          x1s=xx1s[sp]
          for sp2 in [Sparsity.sparse(1,1),Sparsity.dense(1,1),Sparsity.sparse(3,4),Sparsity.dense(3,4), Sparsity(4,3,[0,2,2,3],[1,2,1]).T]:
            for v2 in [0,1,0.2,-0.2]:
              x2_ = DMatrix(sp2,v2)
              xx2 = MX.sym("x",sp2.size1(),sp2.size2())
              x2=xx2[sp2]
              xx2s = SX.sym("x",sp2.size1(),sp2.size2())
              x2s=xx2s[sp2]
              for (casadiop, numpyop,name, flags) in self.matrixbinarypool.zip():
                if "mul" in name and (sp.numel()==1 or sp2.numel()==1): continue
                r = casadiop([x1,x2])
                f = MXFunction([xx1,xx2],[r])
                f.init()
                f.setInput(v1,0)
                f.setInput(v2,1)
                f.evaluate()
                g = MXFunction([xx1,xx2],[r])
                g.init()
                g.setInput(v1,0)
                g.setInput(v2,1)
                g.evaluate()
                
                self.checkarray(f.getOutput(),numpyop([x1_,x2_]),str([sp,sp2,v1,v2,x1_,x2_,name]))
                
                
                if "mul" not in name:
                  a = IMatrix(f.getOutput().sparsity(),1)
                  b = IMatrix(g.getOutput().sparsity(),1)
                  
                  c = b-a
                  if c.size()>0:
                    # At least as sparse as DMatrix calculus
                    self.assertTrue(min(c)>=0,str([sp,sp2,v1,v2,name,a,b]))
                
                if sp.size()>0 and sp2.size()>0 and v1!=0 and v2!=0:
                  self.checkfunction(f,g,hessian=False,failmessage=str([sp,sp2,v1,v2,x1_,x2_,name]))

  @memory_heavy()
  def test_MXConstant(self):
      self.message("MX unary operations, constant")
      
      #self.checkarray(DMatrix(Sparsity.tril(4),1),DMatrix(Sparsity.dense(4,4),1))
      
      for sp in [Sparsity.dense(0,0),Sparsity.dense(0,2),Sparsity.dense(2,0),Sparsity.dense(1,1),Sparsity.dense(2,2), Sparsity(4,3,[0,2,2,3],[1,2,1])]:
        for v in [0,1,0.2]:
          x_ = DMatrix(sp,v)
          
          x=MX(sp,v)
          
          for (casadiop, numpyop,name, flags) in self.pool.zip():
            if 'nozero' in flags and (v==0 or not sp.isDense()): continue
            r = casadiop([x])
            print r
            self.assertTrue(r.isConstant())
            
            self.checkarray(r.getMatrixValue(),numpyop(x_),str([x_,name]))
            
            a = IMatrix(r.getMatrixValue().sparsity(),1)
            b = IMatrix(DMatrix(numpyop(x_)).sparsity(),1)
            
            c = b-a
            if c.size()>0:
              # At least as sparse as DMatrix calculus
              self.assertTrue(min(c)>=0,str([sp,v,name]))
        
      for sp in [Sparsity.dense(1,1),Sparsity.sparse(1,1),Sparsity.sparse(3,4),Sparsity.dense(3,4), Sparsity(4,3,[0,2,2,3],[1,2,1]).T]:
        for v1 in [0,1,0.2,-0.2]:
          x1_ = DMatrix(sp,v1)
          x1=MX(sp,v1)
          for sp2 in [Sparsity.dense(1,1),Sparsity.sparse(1,1),Sparsity.sparse(3,4),Sparsity.dense(3,4), Sparsity(4,3,[0,2,2,3],[1,2,1]).T]:
            for v2 in [0,1,0.2,-0.2]:
              x2_ = DMatrix(sp2,v2)
              x2=MX(sp2,v2)
              for (casadiop, numpyop,name, flags) in self.matrixbinarypool.zip():
                if "mul" in name and (sp.numel()==1 or sp2.numel()==1): continue
                r = casadiop([x1,x2])
                f = MXFunction([],[r]) # Should not be needed -> constant folding
                f.init()
                f.evaluate()
               
                
                self.checkarray(f.getOutput(),numpyop([x1_,x2_]),str([sp,sp2,v1,v2,name]))
                if "mul" not in name:
                  a = IMatrix(f.getOutput().sparsity(),1)
                  b = IMatrix(DMatrix(numpyop([x1_,x2_])).sparsity(),1)
                  
                  c = b-a
                  if c.size()>0:
                    # At least as sparse as DMatrix calculus
                    self.assertTrue(min(c)>=0,str([sp,sp2,v1,v2,name]))

  def test_graph_substitute(self):
    x=MX.sym("X",4,4)
    y=MX.sym("Y",4,4)
    b=MX.sym("B",4,4)

    c = x*y
    d = b*y
    f = c+d
    
    
    C = MX.sym("C",4,4)
    f = graph_substitute(f,[c],[C])
    
    F = MXFunction([y,b,C],[f])
    F.init()
    
    F.setInput(1,0)
    F.setInput(2,1)
    F.setInput(3,2)
    
    F.evaluate()
    
    self.checkarray(F.getOutput(),5*DMatrix.ones(4,4))
    
    D = MX.sym("D",4,4)
    f = graph_substitute(f,[d],[D])
    
    F = MXFunction([D,C],[f])
    F.init()
    
    F.setInput(4,0)
    F.setInput(5,1)
    
    F.evaluate()
    
    self.checkarray(F.getOutput(),9*DMatrix.ones(4,4))
                 

  def test_matrix_expand(self):
    n = 2
    a = MX.sym("a",n,n)
    b = MX.sym("b",n,n)
    c = MX.sym("c",n,n)

    d = a+b
    e = d*c

    self.assertEqual(countNodes(e),6)

    t0 = matrix_expand(e)

    self.assertEqual(countNodes(t0),5)
    
    t1 = matrix_expand(e,[d])
    self.assertEqual(countNodes(t1),6)
    
    print e,t0,t1
    
    
    outs = []
    for x in [e,t0,t1]:
      f = MXFunction([a,b,c],[x])
      f.init()
      
      f.setInput(1.1,0)
      f.setInput(2.2,1)
      f.setInput(3.3,2)
      
      f.evaluate()
      
      outs.append(f.getOutput())
      if outs>1:
        self.checkarray(outs[0],outs[-1])
      
    print outs
    
  def test_kron(self):
    a = sparse(DMatrix([[1,0,6],[2,7,0]]))
    b = sparse(DMatrix([[1,0,0],[2,3,7],[0,0,9],[1,12,13]]))
    
    A = MX.sym("A",a.sparsity())
    B = MX.sym("B",b.sparsity())
    C = c.kron(A,B)
    
    f = MXFunction([A,B],[C])
    f.init()
    f.setInput(a,0)
    f.setInput(b,1)
    f.evaluate()
    
    c_ = f.getOutput()
    
    self.assertEqual(c_.size1(),a.size1()*b.size1())
    self.assertEqual(c_.size2(),a.size2()*b.size2())
    self.assertEqual(c_.size(),a.size()*b.size())
    
    self.checkarray(c_,numpy.kron(a,b))
    
  def test_setSparse(self):
    x = MX.sym("x",Sparsity.tril(3))
    y = x.setSparse(Sparsity.tril(3).T)
    
    f = MXFunction([x],[y])
    f.init()
    
    f.setInput(range(1,4*3/2+1))
    f.evaluate()
    
    self.checkarray(f.getOutput(),DMatrix([[1,0,0],[0,4,0],[0,0,6]]))
    self.checkarray(IMatrix(f.getOutput().sparsity(),1),IMatrix(Sparsity.tril(3).T,1))
    
  def test_repmat(self):
    a = DMatrix([[1,2],[3,4],[5,6]])
    self.checkarray(repmat(a,2,3),kron(DMatrix.ones(2,3),a))
    
  def test_transpose_bug(self):
    
    A = [[-26.9091,00,00,1,00,00,00,00,00,00,00,00,00,00,00],
     [00,-26.9091,00,00,1,00,00,00,00,00,00,00,00,00,00],
     [00,00,-26.9091,00,00,1,00,00,00,00,00,00,00,00,00],
     [00,00,00,-14.1393,-0.348109,-0.382033,3.6491,9.40422,-5.05449,3.03347,-0.00126949,-0.000414104,00,0.00498232,0.00495617],
     [00,00,00,-0.348109,-13.6315,-0.194202,-6.41868,-0.0189287,3.34725,1.86987,-0.000645326,-0.000210504,00,0.0025327,0.0025194],
     [00,00,00,-0.382033,-0.194202,-13.6677,10.1342,-13.0101,-0.890078,-8.99621,-0.000708215,-0.000231018,00,0.00277951,0.00276493],
     [00,00,00,00,00,00,-27.0487,0.36775,-0.277186,0.564995,0.424389,0.0235157,0.227774,00,00],
     [00,00,00,00,00,00,0.113776,-27.3241,0.159254,-0.520988,-0.0235157,0.424389,0.132128,00,00],
     [00,00,00,00,00,00,0.227472,-0.0735535,-26.9135,0.206825,-0.227774,-0.132128,0.424389,00,00],
     [00,00,00,00,00,00,0.332188,-1.02565,-0.0471482,-28.3499,0.132128,-0.227774,0.0235157,00,00],
     [00,00,00,-0.00126949,-0.000645326,-0.000708215,0.00763923,0.0016513,-0.0044564,-0.00284703,-0.119342,-0.000941301,-0.00160961,-1.62792e-05,0.00156863],
     [00,00,00,-0.000414104,-0.000210504,-0.000231018,0.0024919,0.000538651,-0.00145367,-0.000928698,0.000941301,-0.119493,0.000742542,-0.00155303,-7.50989e-06],
     [00,00,00,00,00,00,00,00,00,00,2.1684e-19,-2.1684e-19,-0.171654,-0.000868032,0.000868032],
     [00,00,00,0.00181259,0.000921407,0.0010112,-0.0109074,-0.00235775,0.00636292,0.00406504,0,-0.000567731,-0.000868032,-0.00087529,00],
     [00,00,00,0.00166876,0.000848291,0.000930959,-0.0100419,-0.00217066,0.005858,0.00374247,0.00052268,0,-0.000868032,00,-0.000874062]
     ]
     
    A = sparse(DMatrix(A))

    As = MX.sym("As",A.sparsity())

    f = MXFunction([As],[dense(As.T),dense(As).T,As.T,As,dense(As)])
    f.init()

    f.setInput(A)
    f.evaluate()

    self.checkarray(f.getOutput(0),A.T)
    self.checkarray(f.getOutput(1),A.T)
    self.checkarray(f.getOutput(2),A.T)
    self.checkarray(f.getOutput(3),A)
    self.checkarray(f.getOutput(4),A)
      
  @requiresPlugin(LinearSolver,"csparse")
  def test_bizarre_bug(self):

    A = [[-26.9091,00,00,1,00,00,00,00,00,00,00,00,00,00,00],
     [00,-26.9091,00,00,1,00,00,00,00,00,00,00,00,00,00],
     [00,00,-26.9091,00,00,1,00,00,00,00,00,00,00,00,00],
     [00,00,00,-14.1393,-0.348109,-0.382033,3.6491,9.40422,-5.05449,3.03347,-0.00126949,-0.000414104,00,0.00498232,0.00495617],
     [00,00,00,-0.348109,-13.6315,-0.194202,-6.41868,-0.0189287,3.34725,1.86987,-0.000645326,-0.000210504,00,0.0025327,0.0025194],
     [00,00,00,-0.382033,-0.194202,-13.6677,10.1342,-13.0101,-0.890078,-8.99621,-0.000708215,-0.000231018,00,0.00277951,0.00276493],
     [00,00,00,00,00,00,-27.0487,0.36775,-0.277186,0.564995,0.424389,0.0235157,0.227774,00,00],
     [00,00,00,00,00,00,0.113776,-27.3241,0.159254,-0.520988,-0.0235157,0.424389,0.132128,00,00],
     [00,00,00,00,00,00,0.227472,-0.0735535,-26.9135,0.206825,-0.227774,-0.132128,0.424389,00,00],
     [00,00,00,00,00,00,0.332188,-1.02565,-0.0471482,-28.3499,0.132128,-0.227774,0.0235157,00,00],
     [00,00,00,-0.00126949,-0.000645326,-0.000708215,0.00763923,0.0016513,-0.0044564,-0.00284703,-0.119342,-0.000941301,-0.00160961,-1.62792e-05,0.00156863],
     [00,00,00,-0.000414104,-0.000210504,-0.000231018,0.0024919,0.000538651,-0.00145367,-0.000928698,0.000941301,-0.119493,0.000742542,-0.00155303,-7.50989e-06],
     [00,00,00,00,00,00,00,00,00,00,2.1684e-19,-2.1684e-19,-0.171654,-0.000868032,0.000868032],
     [00,00,00,0.00181259,0.000921407,0.0010112,-0.0109074,-0.00235775,0.00636292,0.00406504,0,-0.000567731,-0.000868032,-0.00087529,00],
     [00,00,00,0.00166876,0.000848291,0.000930959,-0.0100419,-0.00217066,0.005858,0.00374247,0.00052268,0,-0.000868032,00,-0.000874062]
     ]
     
    A = sparse(DMatrix(A))

    b = DMatrix(range(15))
    H = 5

    Bs = MX.sym("Bs",b.sparsity())
    As = MX.sym("As",A.sparsity())

    Ast = As.T

    r= MXFunction([As,Bs],[solve(Ast,Bs,"csparse")])
    r.init()
    R= MXFunction([As,Bs],[solve(dense(Ast),Bs,"csparse")])
    R.init()

    for i in [r,R]:
      i.setInput(A,0)
      i.setInput(b,1)

    r.evaluate()
    R.evaluate()

    r.getOutput().sparsity().spy()
    R.getOutput().sparsity().spy()

    self.checkarray(R.getOutput(),numpy.linalg.solve(A.T,b))
    self.checkarray(r.getOutput(),R.getOutput())

  def test_dependsOn(self):
    a = MX.sym("a")
    b = MX.sym("b")
    
    self.assertTrue(dependsOn(a**2,[a]))
    self.assertTrue(dependsOn(a,[a]))
    self.assertFalse(dependsOn(0,[a]))
    
    self.assertTrue(dependsOn(a**2,[a,b]))
    self.assertTrue(dependsOn(a,[a,b]))
    self.assertFalse(dependsOn(0,[a,b]))
    self.assertTrue(dependsOn(b**2,[a,b]))
    self.assertTrue(dependsOn(b,[a,b]))
    self.assertTrue(dependsOn(a**2+b**2,[a,b]))
    self.assertTrue(dependsOn(a+b,[a,b]))
    self.assertTrue(dependsOn(vertcat([0,a]),[a]))
    self.assertTrue(dependsOn(vertcat([a,0]),[a]))
    self.assertFalse(dependsOn(vertcat([0,b]),[a]))
    self.assertFalse(dependsOn(vertcat([b,0]),[a]))
    self.assertTrue(dependsOn(vertcat([a**2,b**2]),[a,b]))
    self.assertTrue(dependsOn(vertcat([a,0]),[a,b]))
    self.assertTrue(dependsOn(vertcat([0,b]),[a,b]))
    self.assertTrue(dependsOn(vertcat([b,0]),[a,b]))
    self.assertFalse(dependsOn(vertcat([0,0]),[a,b]))
    
  def test_vertcat_simp(self):
    x = MX.sym("x",10)
    y = MX.sym("y")
    z = MX.sym("z")
    x_ = DMatrix(range(10))
    y_ = DMatrix([20])
    z_ = DMatrix([30])
    
    def evalvertcat(a):
      f = MXFunction([x,y,z],[vertcat(a)])
      f.init()
      f.setInput(x_,0)
      f.setInput(y_,1)
      f.setInput(z_,2)
      f.evaluate()
      return f.getOutput()

    self.checkarray(evalvertcat(vertsplit(x)),x_)
    self.checkarray(evalvertcat(vertsplit(x)+[y]),vertcat([x_,y_]))
    self.checkarray(evalvertcat([z]+vertsplit(x)+[y] + vertsplit(x)+[z]),vertcat([z_,x_,y_,x_,z_]))
    self.checkarray(evalvertcat(vertsplit(x)[:-1]),x_[:-1])
    self.checkarray(evalvertcat(vertsplit(x)[:-1]+[y]),vertcat([x_[:-1],y_]))
    self.checkarray(evalvertcat([z]+vertsplit(x)[:-1]+[y] + vertsplit(x)[:-1]+[z]),vertcat([z_,x_[:-1],y_,x_[:-1],z_]))
    self.checkarray(evalvertcat(vertsplit(x)[1:]),x_[1:])
    self.checkarray(evalvertcat(vertsplit(x)[1:]+[y]),vertcat([x_[1:],y_]))
    self.checkarray(evalvertcat([z]+vertsplit(x)[1:]+[y] + vertsplit(x)[1:]+[z]),vertcat([z_,x_[1:],y_,x_[1:],z_]))
    g = vertsplit(x)[5:]+vertsplit(x)[:5]
    self.checkarray(evalvertcat(g),vertcat([x_[5:],x_[:5]]))
    self.checkarray(evalvertcat(g+[y]),vertcat([x_[5:],x_[:5],y_]))
    self.checkarray(evalvertcat([z]+g+[y] + g+[z]),vertcat([z_,x_[5:],x_[:5],y_,x_[5:],x_[:5],z_]))
    
    import __builtin__


    w = vertsplit(x,2)
    r = __builtin__.sum([vertsplit(i) for i in w],[])
    
    self.checkarray(evalvertcat(r),x_)

    w = vertsplit(x,2)
    r = __builtin__.sum([vertsplit(i)+[y] for i in w],[])
    print "vertcat:", r
    print "result:", vertcat(r)

    w = vertsplit(x,2)
    r = __builtin__.sum([vertsplit(i) for i in w],[])
    print "vertcat:", r
    print "result:", vertcat(r+[y])
    
    self.assertTrue(isEqual(vertcat(vertsplit(x)),x))
    
  def test_horzcat_simp(self):
    x = MX.sym("x",1,10)
    y = MX.sym("y")
    z = MX.sym("z")
    x_ = DMatrix(range(10)).T
    y_ = DMatrix([20])
    z_ = DMatrix([30])
    
    def evalhorzcat(a):
      f = MXFunction([x,y,z],[horzcat(a)])
      f.init()
      f.setInput(x_,0)
      f.setInput(y_,1)
      f.setInput(z_,2)
      f.evaluate()
      return f.getOutput()

    self.checkarray(evalhorzcat(horzsplit(x)),x_)
    self.checkarray(evalhorzcat(horzsplit(x)+[y]),horzcat([x_,y_]))
    self.checkarray(evalhorzcat([z]+horzsplit(x)+[y] + horzsplit(x)+[z]),horzcat([z_,x_,y_,x_,z_]))
    self.checkarray(evalhorzcat(horzsplit(x)[:-1]),x_[0,:-1])
    self.checkarray(evalhorzcat(horzsplit(x)[:-1]+[y]),horzcat([x_[0,:-1],y_]))
    self.checkarray(evalhorzcat([z]+horzsplit(x)[:-1]+[y] + horzsplit(x)[:-1]+[z]),horzcat([z_,x_[0,:-1],y_,x_[0,:-1],z_]))
    self.checkarray(evalhorzcat(horzsplit(x)[1:]),x_[0,1:])
    self.checkarray(evalhorzcat(horzsplit(x)[1:]+[y]),horzcat([x_[0,1:],y_]))
    self.checkarray(evalhorzcat([z]+horzsplit(x)[1:]+[y] + horzsplit(x)[1:]+[z]),horzcat([z_,x_[0,1:],y_,x_[0,1:],z_]))
    g = horzsplit(x)[5:]+horzsplit(x)[:5]
    self.checkarray(evalhorzcat(g),horzcat([x_[0,5:],x_[0,:5]]))
    self.checkarray(evalhorzcat(g+[y]),horzcat([x_[0,5:],x_[0,:5],y_]))
    self.checkarray(evalhorzcat([z]+g+[y] + g+[z]),horzcat([z_,x_[0,5:],x_[0,:5],y_,x_[0,5:],x_[0,:5],z_]))
    
    import __builtin__


    w = horzsplit(x,2)
    r = __builtin__.sum([horzsplit(i) for i in w],[])
    
    self.checkarray(evalhorzcat(r),x_)

    w = horzsplit(x,2)
    r = __builtin__.sum([horzsplit(i)+[y] for i in w],[])
    print "vertcat:", r
    print "result:", horzcat(r)

    w = horzsplit(x,2)
    r = __builtin__.sum([horzsplit(i) for i in w],[])
    print "vertcat:", r
    print "result:", horzcat(r+[y])

    self.assertTrue(isEqual(horzcat(horzsplit(x)),x))
    
  def test_vertsplit_simp(self):
    
    dvars = [MX.sym("abcdefghijklm"[i]) for i in range(5) ]
    dvars_ = range(5)

    zz = MX.sym("zz",2)
    zz_ = DMatrix([11,12])
    y = MX.sym("y")
    z = MX.sym("z")
    y_ = DMatrix([20])
    z_ = DMatrix([30])
    
    aa = MX.sym("aa",5)
    aa_ = range(100,105)
    
    def evalvertsplit(a,*args):
      print vertsplit(a,*args)
      f = MXFunction(dvars+[y,z,zz,aa],vertsplit(a,*args))
      f.init()
      for i in range(5):
        f.setInput(dvars_[i],i)
      f.setInput(y_,5+0)
      f.setInput(z_,5+1)
      f.setInput(zz_,5+2)
      f.setInput(aa_,5+3)
      f.evaluate()
      return [f.getOutput(i) for i in range(f.getNumOutputs())]
      
    s= evalvertsplit(vertcat([y]+dvars+[z]))
    self.checkarray(s[0],y_)
    for i in range(5):
      self.checkarray(s[1+i],dvars_[i])
    self.checkarray(s[6],z_)

    s= evalvertsplit(vertcat([y]+dvars+[z]),2)
    
    self.checkarray(s[0],DMatrix([y_,dvars_[0]]))
    self.checkarray(s[1],DMatrix([dvars_[1],dvars_[2]]))
    self.checkarray(s[2],DMatrix([dvars_[3],dvars_[4]]))
    self.checkarray(s[3],DMatrix([z_]))
    
    s= evalvertsplit(vertcat([y,zz,z,zz]),2)
    
    self.checkarray(s[0],DMatrix([y_,zz_[0]]))
    self.checkarray(s[1],DMatrix([zz_[1],z_]))
    self.checkarray(s[2],zz_)
    
    s= evalvertsplit(vertcat([y,zz,z,zz]),3)
    
    self.checkarray(s[0],DMatrix([y_,zz_[0],zz_[1]]))
    self.checkarray(s[1],DMatrix([z_,zz_[0],zz_[1]]))
    
    s= evalvertsplit(vertcat([zz,zz]),2)
    self.checkarray(s[0],zz_)
    self.checkarray(s[1],zz_)

    s= evalvertsplit(vertcat([zz]+dvars))
    self.checkarray(s[0],zz_[0])
    self.checkarray(s[1],zz_[1])
    
    for i in range(5):
      self.checkarray(s[2+i],dvars_[i])

    s= evalvertsplit(vertcat(dvars+[aa]),5)
    self.checkarray(s[0],DMatrix(dvars_))
    self.checkarray(s[1],DMatrix(aa_))

    s= evalvertsplit(vertcat(dvars+[aa]),4)
    self.checkarray(s[0],DMatrix(dvars_[:4]))
    self.checkarray(s[1],DMatrix([dvars_[-1]]+aa_[:3]))
    self.checkarray(s[2],DMatrix(aa_[3:]))

    s= evalvertsplit(vertcat(dvars+[aa]),6)
    self.checkarray(s[0],DMatrix(dvars_+[aa_[0]]))
    self.checkarray(s[1],DMatrix(aa_[1:]))
    
    for i in range(5):
      self.assertTrue(isEqual(vertsplit(vertcat(dvars))[i],dvars[i]))

  def test_horzsplit_simp(self):
    
    dvars = [MX.sym("abcdefghijklm"[i]) for i in range(5) ]
    dvars_ = range(5)

    zz = MX.sym("zz",1,2)
    zz_ = DMatrix([11,12]).T
    y = MX.sym("y")
    z = MX.sym("z")
    y_ = DMatrix([20])
    z_ = DMatrix([30])
    
    aa = MX.sym("aa",1,5)
    aa_ = range(100,105)
    
    def evalhorzsplit(a,*args):
      print horzsplit(a,*args)
      f = MXFunction(dvars+[y,z,zz,aa],horzsplit(a,*args))
      f.init()
      for i in range(5):
        f.setInput(dvars_[i],i)
      f.setInput(y_,5+0)
      f.setInput(z_,5+1)
      f.setInput(zz_,5+2)
      f.setInput(aa_,5+3)
      f.evaluate()
      return [f.getOutput(i) for i in range(f.getNumOutputs())]
      
    s= evalhorzsplit(horzcat([y]+dvars+[z]))
    self.checkarray(s[0],y_)
    for i in range(5):
      self.checkarray(s[1+i],dvars_[i])
    self.checkarray(s[6],z_)

    s= evalhorzsplit(horzcat([y]+dvars+[z]),2)
    
    self.checkarray(s[0],DMatrix([y_,dvars_[0]]).T)
    self.checkarray(s[1],DMatrix([dvars_[1],dvars_[2]]).T)
    self.checkarray(s[2],DMatrix([dvars_[3],dvars_[4]]).T)
    self.checkarray(s[3],DMatrix([z_]).T)
    
    s= evalhorzsplit(horzcat([y,zz,z,zz]),2)
    
    self.checkarray(s[0],DMatrix([y_,zz_[0,0]]).T)
    self.checkarray(s[1],DMatrix([zz_[0,1],z_]).T)
    self.checkarray(s[2],zz_)
    
    s= evalhorzsplit(horzcat([y,zz,z,zz]),3)
    
    self.checkarray(s[0],DMatrix([y_,zz_[0,0],zz_[0,1]]).T)
    self.checkarray(s[1],DMatrix([z_,zz_[0,0],zz_[0,1]]).T)
    
    s= evalhorzsplit(horzcat([zz,zz]),2)
    self.checkarray(s[0],zz_)
    self.checkarray(s[1],zz_)

    s= evalhorzsplit(horzcat([zz]+dvars))
    self.checkarray(s[0],zz_[0,0])
    self.checkarray(s[1],zz_[0,1])
    
    for i in range(5):
      self.checkarray(s[2+i],dvars_[i])

    s= evalhorzsplit(horzcat(dvars+[aa]),5)
    self.checkarray(s[0],DMatrix(dvars_).T)
    self.checkarray(s[1],DMatrix(aa_).T)

    s= evalhorzsplit(horzcat(dvars+[aa]),4)
    self.checkarray(s[0],DMatrix(dvars_[:4]).T)
    self.checkarray(s[1],DMatrix([dvars_[-1]]+aa_[:3]).T)
    self.checkarray(s[2],DMatrix(aa_[3:]).T)

    s= evalhorzsplit(horzcat(dvars+[aa]),6)
    self.checkarray(s[0],DMatrix(dvars_+[aa_[0]]).T)
    self.checkarray(s[1],DMatrix(aa_[1:]).T)
    
    for i in range(5):
      self.assertTrue(isEqual(horzsplit(horzcat(dvars))[i],dvars[i]))
      
  def test_vertsplit_derivative(self):
    m = MX.sym("X",10)

    f = MXFunction([m],[vertsplit(m)[0]])
    f.init()

    f.derivative(0,1)

  def test_MX_const_sp(self):
    x = MX.sym("x",4,1)

    sp = Sparsity.triplet(3,3,[0,1,2,2],[0,0,1,2])

    f = MXFunction([x],[x[IMatrix(sp,range(sp.size()))]])
    f.init()

    g = MXFunction([x],[MX(sp,x)])
    g.init()
    
    f.setInput(range(1,5))
    g.setInput(range(1,5))
    
    self.checkfunction(f,g)
    
  def test_reshape_sp(self):
    x = MX.sym("x",4,1)

    f = MXFunction([x],[x.reshape((2,2))])
    f.init()
    
    sx = SX.sym("x",4,1)

    g = SXFunction([sx],[sx.reshape((2,2))])
    g.init()
    
    f.setInput(range(1,5))
    g.setInput(range(1,5))
    
    self.checkfunction(f,g)
    
  def test_issue1041(self):
    x = MX.sym("x",2)

    y = vertsplit(x,[0,1,2])[1]

    f = MXFunction([x],[y])
    f.init()

    H = f.hessian()
    H.init()
    
  def test_bug_1042(self):

    x = MX.sym('x',2,1)
    
    mf = MXFunction([x],[x*x[0,0]])
    mf.init()
    
    mfunction = mf.expand()
    mfunction.init()
    
    mfg = mf.derivative(0,1)
    mfg.init()
    
    mfunctiong = mfunction.derivative(0,1)
    mfunctiong.init()
    
    for f in [mfg,mfunctiong]:
      f.setInput([1,2],0)
      #f.setInput(0.1,)
      f.setInput([4,5],1)
    
    self.checkfunction(mfg,mfunctiong)
    
  def test_bug_1042bis(self):
    x = MX.sym('x',2,1)
    a = MX.sym("ax",2,1)
    i1 = x[0,0]
    z = i1*x
    i3 = i1*a
    i3= inner_prod(x,a)
    d = MXFunction([x,a],[z,i3])
    d.init()
    d.setInput([1,2],0)
    d.setInput([3,4],1)

    dx = d.expand()
    dx.init()
    dx.setInput([1,2],0)
    dx.setInput([3,4],1)
    
    self.checkfunction(d,dx)
    
  def test_bug_1042tris(self):
    x = MX.sym('x',2,1)
    a = MX.sym("ax",2,1)
    d = MXFunction([x,a],[inner_prod(x,a)])
    d.init()
    d.setInput([1,2],0)
    d.setInput([3,4],1)

    dx = d.expand()
    dx.init()
    dx.setInput([1,2],0)
    dx.setInput([3,4],1)
    
    self.checkfunction(d,dx)
    
  def test_bug_1046(self):
    x = MX.sym('x',1,1)
    y = MX.sym('y',1,1)
    z = jacobian(x,y)
    
    self.assertTrue(z.size()==0)
    
  def test_singularcat(self):

    for c in [MX,SX,DMatrix]:
      x0 = c.zeros(10,0)
      
      x1s = vertsplit(x0, [0,5,10])
      
      for x in x1s:
        self.checkarray(x.shape,(5,0))


      x2 = vertcat(x1s)
      self.checkarray(x2.shape,(10,0))
      
      x2 = vertcat([c.zeros(0,0)] + x1s + [c.zeros(0,0)])
      self.checkarray(x2.shape,(10,0))
        
    for c in [MX,SX,DMatrix]:
      x0 = c.zeros(0,10)
      
      x1s = horzsplit(x0, [0,5,10])
      
      for x in x1s:
        self.checkarray(x.shape,(0,5))

      x2 = horzcat(x1s)
      self.checkarray(x2.shape,(0,10))
      
      x2 = horzcat([c.zeros(0,0)] + x1s + [c.zeros(0,0)])
      self.checkarray(x2.shape,(0,10))
 
    for c in [MX,SX,DMatrix]:
      x0 = c.zeros(10,0)
      
      x1s = vertsplit(x0, [0,5,10])

      x0 = c.zeros(0,10)      
      x1st = horzsplit(x0, [0,5,10])
      
      x2 = blkdiag(x1s)
      self.checkarray(x2.shape,(10,0))
      
      x2 = blkdiag([c.zeros(0,0)] + x1s + [c.zeros(0,0)])
      self.checkarray(x2.shape,(10,0))

      x2 = blkdiag(x1st)
      self.checkarray(x2.shape,(0,10))
      
      x2 = blkdiag([c.zeros(0,0)] + x1st + [c.zeros(0,0)])
      self.checkarray(x2.shape,(0,10))
      
      x2 = blkdiag(x1s+x1st)
      self.checkarray(x2.shape,(10,10))
      
      x2 = blkdiag([c.zeros(0,0)] + x1s+x1st + [c.zeros(0,0)])
      self.checkarray(x2.shape,(10,10))
  def test_empty_symm_jac(self):

    x = MX.sym("x",2)

    g = MXFunction([x],[MX.sparse(1,1)])
    g.init()

    h = g.jacobian(0,0,False,True)
      
    x = MX.sym("x",2)

    g = MXFunction([x],[MX.zeros(1,1)])
    g.init()

    h = g.jacobian(0,0,False,True)

if __name__ == '__main__':
    unittest.main()
