from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *
from scipy.sparse import *

class Matrixtests(casadiTestCase):
  def setUp(self):
    """
    We use this sparse matrix as the test object:
    | 1 0 0 0 0 0 0 |
    | 1 0 0 1 0 0 0 |
    | 0 0 1 0 0 0 0 |
    | 0 0 0 0 0 0 1 |
    | 0 0 0 0 0 0 0 |
    | 1 1 1 1 1 1 1 |
    """
    self.nz = [(0,0),(1,0),(1,3),(2,2),(3,6),(5,0),(5,1),(5,2),(5,3),(5,4),(5,5),(5,6)]
    self.A = DMatrix(6,7)
    self.a = lil_matrix((6,7))
    for p in self.nz:
      self.A[p[0],p[1]]=1
      self.a[p[0],p[1]]=1
      
      
  def test_size12(self):
    self.assertEqual(self.A.size1(),self.a.shape[0],"size")
    self.assertEqual(self.A.size2(),self.a.shape[1],"size")
    
  def test_numel(self):
    self.assertEqual(self.A.numel(),self.A.size1()*self.A.size2(),"numel")

  def test_size(self):
    self.assertEqual(self.A.size(),self.a.nnz,"size")

  def test_sizeU(self):
    a = 0
    for p in self.nz:
      if p[1]>=p[0]:
        a=a+1
    self.assertEqual(self.A.sparsity().sizeU(),a,"sizeU")
    
  def test_sizeL(self):
    a = 0
    for p in self.nz:
      if p[1]<=p[0]:
        a=a+1
    self.assertEqual(self.A.sparsity().sizeL(),a,"sizeL")
    
  def test_col(self):
    col=list(self.A.sparsity().col())
    for k in range(self.A.size()):
      self.assertEqual(col[k],self.nz[k][1],"col")

  def test_row0(self):
    row=list(self.A.sparsity().getRow())
    for k in range(self.A.size()):
      self.assertEqual(row[k],self.nz[k][0],"row")
      
  def test_row(self):
    row=list(self.A.sparsity().rowind())
    for k in range(self.A.size()):
      self.assertTrue(row[self.nz[k][0]] <= k,"row")
      self.assertTrue(k<=row[self.nz[k][0]+1],"row")
      
  def test_col2(self):
    for k in range(self.A.size()):
      self.assertEqual(self.A.sparsity().col(k),self.nz[k][1],"col2")
   
  def test_row2(self):
    for k in range(self.A.size()):
      self.assertTrue(self.A.sparsity().rowind(self.nz[k][0]) <= k,"row2")
      self.assertTrue(k<=self.A.sparsity().rowind(self.nz[k][0]+1),"row2")

  def test_getNZ(self):
    k=0
    for i in range(self.A.size1()):
      for j in range(self.A.size2()):
        if (i,j) in self.nz:
          self.assertEquals(self.A.sparsity().getNZ(i,j),k,"getNZ")
          k=k+1
        else:
          # how to force the getN const to be called?
          #self.assertEquals(self.A.sparsity().getNZ(i,j),-1,"getNZ")
          pass
   
if __name__ == '__main__':
    unittest.main()
