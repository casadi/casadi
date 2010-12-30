from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *

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
    zt = f.getOutput()
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
    zt = f.getOutput()
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
    zt = f.getOutput()
    zr = array([[L[0]*i,L[1]*i,L[2]*i] for i in range(8)])
    checkarray(self,zrf(zr),zt,name)
    return (zt,zrf(zr))
    
class MXtests(unittest.TestCase):

  def setUp(self):
    pass

  def test_MX1(self):
    x = MX("x",2,3)
    self.assertEqual(x.size1(),2,"MX fails to indicate its size1")
    self.assertEqual(x.size2(),3,"MX fails to indicate its size2")

  def test_MXvertcat(self):
    x = MX("x",1,3)
    y = MX("y",1,3)
    z=vertcat(x,y)
    self.assertEqual(z.size1(),2,"MX fails to indicate its size1")
    self.assertEqual(z.size2(),3,"MX fails to indicate its size2")

  def test_MXindices(self):
    x = MX("x",2,3)
    for i in range(2):
      for j in range(3):
        a = "x(%d,%d)" % (i,j)
        b = str(x[i,j])
        self.assertEqual(a,b,"MX indexing is mixed up")

  def test_MXslicing(self):
    x = MX("x",2,3)
    y=x.getColumn(0)
    self.assertEqual(y.size1(),2,"GetColumn returns MX of wrong dimension")
    self.assertEqual(y.size2(),1,"GetColumn returns MX of wrong dimension")

    self.assertEqual(str(y[0]),"[[x(0,0),];[x(1,0),];][0]","Slicing error")
    self.assertEqual(str(y[1]),"[[x(0,0),];[x(1,0),];][1]","Slicing error")

    z=x.getRow(1)
    self.assertEqual(z.size1(),1,"GetRow returns MX of wrong dimension")
    self.assertEqual(z.size2(),3,"GetRow returns MX of wrong dimension")

    self.assertEqual(str(z[0]),"[[x(1,0),x(1,1),x(1,2),];](0,0)","Slicing error")
    self.assertEqual(str(z[1]),"[[x(1,0),x(1,1),x(1,2),];](0,1)","Slicing error")
    self.assertEqual(str(z[2]),"[[x(1,0),x(1,1),x(1,2),];](0,2)","Slicing error")

    x.getRow(0)
    x.getColumn(0)
    x.getColumn(1)
    x.getColumn(2)

  def test_MXslicing2(self):
    x = MX("x",2,3)
    y=x[:,0]
    self.assertEqual(y.size1(),2,"GetColumn returns MX of wrong dimension")
    self.assertEqual(y.size2(),1,"GetColumn returns MX of wrong dimension")

    self.assertEqual(str(y[0]),"[[x(0,0),];[x(1,0),];][0]","Slicing error")
    self.assertEqual(str(y[1]),"[[x(0,0),];[x(1,0),];][1]","Slicing error")

    z=x[1,:]
    self.assertEqual(z.size1(),1,"GetRow returns MX of wrong dimension")
    self.assertEqual(z.size2(),3,"GetRow returns MX of wrong dimension")

    self.assertEqual(str(z[0]),"[[x(1,0),x(1,1),x(1,2),];](0,0)","Slicing error")
    self.assertEqual(str(z[1]),"[[x(1,0),x(1,1),x(1,2),];](0,1)","Slicing error")
    self.assertEqual(str(z[2]),"[[x(1,0),x(1,1),x(1,2),];](0,2)","Slicing error")

  def test_MXslicing3(self):
    x = MX("x",2,3)
    y=x[[0,1],0]
    self.assertEqual(y.size1(),2,"GetColumn returns MX of wrong dimension")
    self.assertEqual(y.size2(),1,"GetColumn returns MX of wrong dimension")

    self.assertEqual(str(y[0]),"[[x(0,0),];[x(1,0),];][0]","Slicing error")
    self.assertEqual(str(y[1]),"[[x(0,0),];[x(1,0),];][1]","Slicing error")

    z=x[1,[0,1,2]]
    self.assertEqual(z.size1(),1,"GetRow returns MX of wrong dimension")
    self.assertEqual(z.size2(),3,"GetRow returns MX of wrong dimension")

    self.assertEqual(str(z[0]),"[[x(1,0),x(1,1),x(1,2),];](0,0)","Slicing error")
    self.assertEqual(str(z[1]),"[[x(1,0),x(1,1),x(1,2),];](0,1)","Slicing error")
    self.assertEqual(str(z[2]),"[[x(1,0),x(1,1),x(1,2),];](0,2)","Slicing error")


  def test_MXslicing4(self):
    x = MX("x",4,5)
    y=x[[0,1],[1,4,3]]
    self.assertEqual(y.size1(),2,"GetColumn returns MX of wrong dimension")
    self.assertEqual(y.size2(),3,"GetColumn returns MX of wrong dimension")

  def test_MXslicing5(self):
    x = MX("x",2,3)
    z=x[1,0:-1]
    self.assertEqual(z.size1(),1,"GetRow returns MX of wrong dimension")
    self.assertEqual(z.size2(),3,"GetRow returns MX of wrong dimension")

    self.assertEqual(str(z[0]),"[[x(1,0),x(1,1),x(1,2),];](0,0)","Slicing error")
    self.assertEqual(str(z[1]),"[[x(1,0),x(1,1),x(1,2),];](0,1)","Slicing error")
    self.assertEqual(str(z[2]),"[[x(1,0),x(1,1),x(1,2),];](0,2)","Slicing error")

  def test_MXslicing6(self):
    x = MX("x",2,3)
    z=x[1,0:-2]
    self.assertEqual(z.size1(),1,"GetRow returns MX of wrong dimension")
    self.assertEqual(z.size2(),2,"GetRow returns MX of wrong dimension")
    self.assertEqual(str(z[0]),"[[x(1,0),x(1,1),];](0,0)","Slicing error")
    self.assertEqual(str(z[1]),"[[x(1,0),x(1,1),];](0,1)","Slicing error")

  def test_MXslicing7(self):
    x = MX("x",2,3)
    z=x[1,::2]
    self.assertEqual(z.size1(),1,"GetRow returns MX of wrong dimension")
    self.assertEqual(z.size2(),2,"GetRow returns MX of wrong dimension")
    self.assertEqual(str(z[0]),"[[x(1,0),x(1,2),];](0,0)","Slicing error")
    self.assertEqual(str(z[1]),"[[x(1,0),x(1,2),];](0,1)","Slicing error")

  def test_MX2(self):
    U = MX("U",10,2)
    u = U.getRow(1)

  def test_MXfunction1(self):
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
    yt = f.getOutputData()
    self.assertEqual(type(yt),TupleType,"Output of MXFunction is expected to be tuple of floats")
    self.assertEqual(len(yt),1,"Output of MXFunction was tuple of floats, as expected, but length is incorrect.")
    y=yt[0]
    self.assertEqual(type(y),float,"Output of MXFunction is expected to be tuple of floats")
    self.assertAlmostEqual(y, 2*3,10)

  def test_MXfunction2(self):
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
    zt1 = f.getOutputData(0)
    zt2 = f.getOutputData(1)
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
    zt1 = f.getOutputData(0)
    zt2 = f.getOutputData(1)
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
    zt1 = f.getOutput(0)
    zt2 = f.getOutput(1)
    
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
    zt=f.getOutput()
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
    zt = f.getOutput()
    self.assertEqual(type(zt),ndarray,"Output of MXFunction is expected to be numpy.ndarray")
    self.assertEqual(zt.shape[0],1,"Output of MXFunction is of wrong shape.")
    self.assertEqual(zt.shape[1],2,"Output of MXFunction is of wrong shape.")
    z1=zt[0,0]
    z2=zt[0,1]
    self.assertEqual(type(z1),float64,"Output of MXFunction is expected to be numpy.ndarray of floats")
    self.assertEqual(type(z2),float64,"Output of MXFunction is expected to be numpy.ndarray of floats")
    self.assertAlmostEqual(z2, 21,10)
    self.assertAlmostEqual(z1, 10,10)

  def test_MXslice(self):
    x = MX("x",1,3)
    z=x[:]
    self.assertEqual(z.size1(),1,"Flatten returns MX of wrong dimension")
    self.assertEqual(z.size2(),3,"Flatten returns MX of wrong dimension")
    f = MXFunction([x],[z])
    f.init()
    L=[1,2,3]
    f.setInput(L,0)
    f.evaluate()
    zt = f.getOutput()
    
    for i in range(3):
      self.assertAlmostEqual(L[i], zt[i],10)
      
  def test_MXslice(self):
    x = MX("x",3,1)
    z=x[:]
    self.assertEqual(z.size1(),3,"Flatten returns MX of wrong dimension")
    self.assertEqual(z.size2(),1,"Flatten returns MX of wrong dimension")
    f = MXFunction([x],[z])
    f.init()
    L=[1,2,3]
    f.setInput(L,0)
    f.evaluate()
    zt = f.getOutput()
    
    for i in range(3):
      self.assertAlmostEqual(L[i], zt[i],10)
        
  def test_MXorder(self):
    x = MX("x",2,3)
    f = MXFunction([x],[x])
    self.assertEqual(f.getNumInputs(),1,"MXFunction fails to indicate correct number of inputs")
    self.assertEqual(f.getNumOutputs(),1,"MXFunction fails to indicate correct number of outputs")
    f.init()
    L=[1,2,3,4,5,6]
    f.setInput(L,0)
    f.evaluate()
    zt = f.getOutput()
    self.assertEqual(zt.shape[0],2,"Output of MXFunction is of wrong shape.")
    self.assertEqual(zt.shape[1],3,"Output of MXFunction is of wrong shape.")
      
    Lr=reshape(L,(2,3))
    for i in range(2):
      for j in range(3):
        self.assertAlmostEqual(Lr[i,j], zt[i,j],10)
    
  def test_MXtrans(self):
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
    zt = f.getOutput()
    
    ztr=reshape(zt,(3,2))
    Lr=reshape(L,(2,3))
    for i in range(2):
      for j in range(3):
        self.assertAlmostEqual(Lr[i,j], ztr[j,i],10)
    
  def test_MXflatten(self):
    x = MX("x",2,3)
    z=flatten(x)
    self.assertEqual(z.size1(),6,"Flatten returns MX of wrong dimension")
    self.assertEqual(z.size2(),1,"Flatten returns MX of wrong dimension")
    f = MXFunction([x],[z])
    self.assertEqual(f.getNumInputs(),1,"MXFunction fails to indicate correct number of inputs")
    self.assertEqual(f.getNumOutputs(),1,"MXFunction fails to indicate correct number of outputs")
    f.init()
    L=[1,2,3,4,5,6]
    f.setInput(L,0)
    f.evaluate()
    zt = f.getOutput()
    for i in range(len(L)):
      self.assertAlmostEqual(L[i], zt[i],10)
      
    
  def test_MXreshape(self):
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
    zt = f.getOutput()
    for i in range(len(L)):
      self.assertAlmostEqual(L[i], zt[0,i],10)
  
  def test_MXcompose(self):
    checkMXoperations(self,lambda x: x,lambda x: x,'vertcat')
    checkMXoperations(self,lambda x: trans(x),lambda x: x.T,'trans(vertcat)')
    checkMXoperations(self,lambda x: trans(trans(x)),lambda x: x,'trans(trans(vertcat))')
    checkMXoperations(self,lambda x: flatten(trans(x)),lambda x: reshape(x.T,(prod(x.shape),1)),'flatten(trans(vertcat))')
    checkMXoperations(self,lambda x: trans(flatten(x)),lambda x: reshape(x,(prod(x.shape),1)).T,'flatten(trans(vertcat))')
    checkMXoperations(self,lambda x: c.reshape(x,(4,6)),lambda x: reshape(x,(4,6)),'reshape(vertcat)')
    checkMXoperations(self,lambda x: c.reshape(trans(x),(4,6)),lambda x: reshape(x.T,(4,6)),'reshape(trans(vertcat))') 
    checkMXoperations(self,lambda x: trans(c.reshape(x,(4,6))),lambda x: reshape(x,(4,6)).T,'trans(reshape(vertcat))') 

  def test_MXcompose2(self):
    checkMXoperations2(self,lambda x: x,lambda x: x,'horzcat')
    checkMXoperations2(self,lambda x: trans(x),lambda x: x.T,'trans(horzcat)')
    checkMXoperations2(self,lambda x: trans(trans(x)),lambda x: x,'trans(trans(horzcat))')
    checkMXoperations2(self,lambda x: flatten(trans(x)),lambda x: reshape(x.T,(prod(x.shape),1)),'flatten(trans(horzcat))')
    checkMXoperations2(self,lambda x: trans(flatten(x)),lambda x: reshape(x,(prod(x.shape),1)).T,'flatten(trans(horzcat))')
    checkMXoperations2(self,lambda x: c.reshape(x,(4,6)),lambda x: reshape(x,(4,6)),'reshape(horzcat)')
    checkMXoperations2(self,lambda x: c.reshape(trans(x),(4,6)),lambda x: reshape(x.T,(4,6)),'reshape(trans(horzcat))') 
    checkMXoperations2(self,lambda x: trans(c.reshape(x,(4,6))),lambda x: reshape(x,(4,6)).T,'trans(reshape(horzcat))') 

  def test_MXcompose3(self):
    checkMXoperations3(self,lambda x: x,lambda x: x,'snippet')
    checkMXoperations3(self,lambda x: trans(x),lambda x: x.T,'trans(snippet)')
    checkMXoperations3(self,lambda x: trans(trans(x)),lambda x: x,'trans(trans(snippet))')
    checkMXoperations3(self,lambda x: flatten(trans(x)),lambda x: reshape(x.T,(prod(x.shape),1)),'flatten(trans(snippet))')
    checkMXoperations3(self,lambda x: trans(flatten(x)),lambda x: reshape(x,(prod(x.shape),1)).T,'flatten(trans(snippet))')
    checkMXoperations3(self,lambda x: c.reshape(x,(4,6)),lambda x: reshape(x,(4,6)),'reshape(snippet)')
    checkMXoperations3(self,lambda x: c.reshape(trans(x),(4,6)),lambda x: reshape(x.T,(4,6)),'reshape(trans(snippet))') 
    checkMXoperations3(self,lambda x: trans(c.reshape(x,(4,6))),lambda x: reshape(x,(4,6)).T,'trans(reshape(snippet))') 

  def test_MXcompose4(self):
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
    x = MX("x",2,3)
    xt=trans(x)
    self.assertEqual(xt.size1(),3,"trans MX of wrong dimension")
    self.assertEqual(xt.size2(),2,"trans MX of wrong dimension")
    p1 = horzcat([xt[0,0],xt[1,0],xt[2,0]])
    p2 = horzcat([xt[0,1],xt[1,1],xt[2,1]])
    z = vertcat([p1,p2])
    f = MXFunction([x],[z])
    f.init()
    L=[1,2,3,4,5,6]
    f.setInput(L,0)
    f.evaluate()
    zt = f.getOutput()
    zr = reshape(array(L),(2,3))
    checkarray(self,zr,zt,"slicing(trans)")
    checkMXoperations(self,lambda x: x[:,0],lambda x: x[:,0],'vertcat[:,0]')
    checkMXoperations(self,lambda x: x[:,1],lambda x: x[:,1],'vertcat[:,1]')
    checkMXoperations(self,lambda x: x[1,:],lambda x: x[1,:],'vertcat[1,:]')
    checkMXoperations(self,lambda x: x[0,:],lambda x: x[0,:],'vertcat[0,:]')
    checkMXoperations(self,lambda x: x[:,0:1],lambda x: x[:,0:1],'vertcat[:,0:1]')
    checkMXoperations(self,lambda x: x[0:1,:],lambda x: x[0:1,:],'vertcat[0:1,:]')
    checkMXoperations(self,lambda x: x[[0,1],0:1],lambda x: x[[0,1],0:1],'vertcat[:,0:1]')
    checkMXoperations(self,lambda x: x[0:1,[0,1]],lambda x: x[0:1,[0,1]],'vertcat[0:1,:]')
    
if __name__ == '__main__':
    unittest.main()

