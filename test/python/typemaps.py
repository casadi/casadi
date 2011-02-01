from casadi import *
import casadi as c
from numpy import *
import numpy as n
import unittest
from types import *
from helpers import *
from scipy.sparse import *
from scipy import *
import warnings

class typemaptests(casadiTestCase):

  def setUp(self):
    pass

  def test_0(self):
    arrays = [array([[1,2,3],[4,5,6]]),array([[1,2],[3,4],[5,6]],dtype=double),array([[3.2,4.6,9.9]])]
    for i in range(len(arrays)):
      m = arrays[i]
      zt=trans(trans(m))
      self.assertTrue(isinstance(zt,DMatrix),"DMatrix expected")
      self.checkarray(m,zt,"DMatrix(numpy.ndarray)")
      self.checkarray(m,zt.toArray(),"DMatrix(numpy.ndarray).toArray()")
      self.checkarray(m,zt.toCsr_matrix(),"DMatrix(numpy.ndarray).toCsr_matrix()")
      
  def test_1(self):
    arrays = [DMatrix(3,4,[1,2,1],[0,2,2,3],[3,2.3,8])]
    for i in range(len(arrays)):
      m = arrays[i]
      zt=trans(trans(m))
      self.assertTrue(isinstance(zt,DMatrix),"DMatrix expected")
      self.checkarray(m,zt,"DMatrix(DMatrix)")
      self.checkarray(m,zt.toArray(),"DMatrix(DMatrix).toArray()")
      self.checkarray(m,zt.toCsr_matrix(),"DMatrix(DMatrix).toCsr_matrix()")
   
  def test_2(self):
    arrays = [csr_matrix( ([3,2.3,8],([0,2,0],[1,1,2])), shape = (3,4), dtype=double ),
              csr_matrix( ([3,2.3,8],([0,2,0],[1,1,2])), shape = (3,4), dtype=int )
              ]
    for i in range(len(arrays)):
      m = arrays[i]
      zt=trans(trans(m))
      self.assertTrue(isinstance(zt,DMatrix),"DMatrix expected")
      self.checkarray(m,zt,"DMatrix(crs_matrix)")
      self.checkarray(m,zt.toArray(),"DMatrix(crs_matrix).toArray()")
      self.checkarray(m,zt.toCsr_matrix(),"DMatrix(crs_matrix).toCsr_matrix()")
      
      
  def test_setget(self):
    data = n.array([3,2.3,8])
    dm=DMatrix(3,4,[1,2,1],[0,2,2,3],[3,2.3,8])
    
    c=dm.toCsr_matrix()
    z=n.zeros((3,4))
    dm.get(z)
    self.checkarray(z,dm,"get(2Dndarray)")
    z=n.matrix(n.zeros((3,4)))
    dm.get(z)
    self.checkarray(z,dm,"get(2Dmatrix)")
    z=n.zeros((12,5))
    self.assertRaises(TypeError,lambda : dm.get(z),"get(wrong size ndarray)")
    z=ones((3,4))
    dm.set(z)
    self.checkarray(dm.toArray() > 0,dm,"set(2Dndarray)")
    z=n.matrix(ones((3,4)))
    dm.set(z)
    self.checkarray(dm.toArray() > 0,dm,"set(2Dmatrix)")
    z=n.zeros((12,5))
    self.assertRaises(TypeError,lambda : dm.set(z))
    
    dm.set(c)
    self.checkarray(c,dm,"set(csr_matrix)")
    
    z=n.zeros(3)
    dm.get(z)
    self.checkarray(n.matrix(z),n.matrix(data),"get(1Dndarray)")
    dm.set(z)
    self.checkarray(c,dm,"set(1Dndarray)")

    dm = dm * 2
    dm.get(c)
    dm.shape = (dm.size1(),dm.size2())

    self.checkarray(c,dm,"get(csr_matrix)")
    
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      c[0,0]=1
    self.assertRaises(TypeError,lambda :  dm.set(c))
    self.assertRaises(TypeError,lambda :  dm.get(c))

  def test_conversion(self):
    w = DMatrix(3,4,[1,2,1],[0,2,2,3],[3,2.3,8])
    d = array([[1,2,3],[4,5,6]])
    
    list(w)
    tuple(w)
    w.toArray()
    array(w)
    w.toMatrix()
    matrix(w)
    w.toCsr_matrix()

    self.checkarray(DMatrix(d),d,"DMatrix(numpy.ndarray)")
    #self.checkarray(DMatrix(array([1,2,3,4,5,6])),d.ravel(),"DMatrix(numpy.ndarray)")
    #print DMatrix(array([1,2,3,6]),2,2).toArray()

    #print DMatrix(array([1,2,3,6]),2,2).toArray()

if __name__ == '__main__':
    unittest.main()
