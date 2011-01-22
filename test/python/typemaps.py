from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *
from scipy.sparse import *
from scipy import *
    
class typemaptests(casadiTestCase):

  def setUp(self):
    pass

  def test_0(self):
    arrays = [array([[1,2,3],[4,5,6]]),array([[1,2],[3,4],[5,6]],dtype=double),array([[3.2,4.6,9.9]])]
    for i in range(len(arrays)):
      m = arrays[i]
      zt=trans(trans(m))
      self.assertTrue(isinstance(zt,DMatrix),"DMatrix expected")
      zt.shape = (zt.size1(),zt.size2())
      self.checkarray(m,zt,"DMatrix(numpy.ndarray)")
      self.checkarray(m,zt.toArray(),"DMatrix(numpy.ndarray).toArray()")
      self.checkarray(m,zt.toCsr_matrix(),"DMatrix(numpy.ndarray).toCsr_matrix()")
      
  def test_1(self):
    arrays = [DMatrix(3,4,[1,2,1],[0,2,2,3],[3,2.3,8])]
    for i in range(len(arrays)):
      m = arrays[i]
      zt=trans(trans(m))
      self.assertTrue(isinstance(zt,DMatrix),"DMatrix expected")
      zt.shape = (zt.size1(),zt.size2())
      m.shape = (m.size1(),m.size2())
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
      zt.shape = (zt.size1(),zt.size2())
      self.checkarray(m,zt,"DMatrix(crs_matrix)")
      self.checkarray(m,zt.toArray(),"DMatrix(crs_matrix).toArray()")
      self.checkarray(m,zt.toCsr_matrix(),"DMatrix(crs_matrix).toCsr_matrix()")
      
      
  def test_setget(self):
    dm=DMatrix(3,4,[1,2,1],[0,2,2,3],[3,2.3,8])
    dm.shape = (dm.size1(),dm.size2())
    
    c=dm.toCsr_matrix()
    
    z=zeros((3,4))
    dm.get(z)
    self.checkarray(z,dm,"get(ndarray)")
    z=ones((3,4))
    dm.set(z)
    self.checkarray(dm.toArray() > 0,dm,"set(ndarray)")
    
    dm.set(c)
    self.checkarray(c,dm,"set(csr_matrix)")
    
    dm = dm * 2
    dm.get(c)
    dm.shape = (dm.size1(),dm.size2())

    self.checkarray(c,dm,"get(csr_matrix)")
      
    z=zeros((4,3))
    #self.assertRaises(lambda : dm.get(z), Exception)
    
    z=list(zeros(12))
    #print z
    #dm.get(z)
    #print z
    
if __name__ == '__main__':
    unittest.main()

