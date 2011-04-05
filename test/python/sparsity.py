from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *

class Sparsitytests(casadiTestCase):
  def test_union(self):
    self.message("Sparsity union")
    nza = set([  (0,0),
             (0,1),
             (2,0),
             (3,1)])
    nzb = set([  (0,2),
             (0,0),
             (2,2)])
    
    a = CRSSparsity(4,5)
    for i in nza:
      a.getNZ(i[0],i[1])
      
    b = CRSSparsity(4,5)  
    for i in nzb:
      b.getNZ(i[0],i[1])
      
    w = IVector()
    c=a.patternUnion(b,w)
    self.assertEquals(w.size(),len(nza.union(nzb)))
    for k in range(w.size()):
      ind = (c.getRow()[k],c.col(k))
      if (ind in nza and ind in nzb):
        self.assertEquals(w[k],0)
      elif (ind in nza):
        self.assertEquals(w[k],-1)
      elif (ind in nzb):
        self.assertEquals(w[k],1)
        
  def test_getNZDense(self):
    self.message("getNZDense")
    nza = set([  (0,0),(0,1),(2,0),(3,1)])
    
    a = CRSSparsity(4,5)
    for i in nza:
      a.getNZ(i[0],i[1])
      
    A = DMatrix(a,1)
    Ad = DMatrix(array(A))
    for i in getNZDense(a):
      self.assertEqual(Ad[i],1)

  def test_enlarge(self):
    self.message("enlarge")
    import numpy
    self.message(":dense")
    #sp = CRSSparsity(3,4,[1,2,1],[0,2,2,3])
    sp = CRSSparsity(3,4,True)
    
    col = [1,2,4]
    row = [0,3,4,6]
    sp.enlarge(7,8,col,row)
    
    z = numpy.zeros((7,8))
    for i in col:
      for j in row:
        z[i,j]=1

    self.checkarray(DMatrix(sp,1),z,"enlarge")
    self.message(":sparse")
    sp = CRSSparsity(3,4,[1,2,1],[0,2,2,3])
    n = DMatrix(sp,1)
    z = numpy.zeros((7,8))
    for i in range(3):
      for j in range(4):
          z[col[i],row[j]]= n[i,j]
    sp.enlarge(7,8,[1,2,4],[0,3,4,6])
    
    self.checkarray(DMatrix(sp,1),z,"enlarge")
    
if __name__ == '__main__':
    unittest.main()

