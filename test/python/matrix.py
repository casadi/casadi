from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *
import numpy

class Matrixtests(casadiTestCase):
  def test_constructorlol(self):
    self.message("List of list constructor")
    a=DMatrix(array([[1,2,3],[4,5,6],[7,8,9]]))
    b=DMatrix([[1,2,3],[4,5,6],[7,8,9]])
    self.checkarray(a,b,"List of list constructor")
    
  def test_sum(self):
    self.message("sum")
    D=DMatrix([[1,2,3],[4,5,6],[7,8,9]])
    self.checkarray(c.sumRows(D),array([[12,15,18]]),'sum()')
    self.checkarray(c.sumCols(D),array([[6,15,24]]).T,'sum()')
    
  def test_inv(self):
    self.message("Matrix inverse")
    a = DMatrix([[1,2],[1,3]])
    self.checkarray(mul(c.inv(a),a),eye(2),"DMatrix inverse")

  def test_iter(self):
    self.message("iterator")
    L = []
    for i in DMatrix([5,6,7,8]):
      L.append(i)
    self.assertEquals(L[0],5)
    self.assertEquals(L[1],6)
    self.assertEquals(L[2],7)
    self.assertEquals(L[3],8)
    
  def test_tuple_unpacking(self):
    self.message("tuple unpacking")
    (a,b,c,d) = DMatrix([5,6,7,8])
    self.assertEquals(a,5)
    self.assertEquals(b,6)
    self.assertEquals(c,7)
    self.assertEquals(d,8)
    
  def test_numpy(self):
    self.message("numpy check")
    # This is an example that failed on a windows machine
    import numpy as NP
    A = NP.zeros((3,4),dtype=SX)

    x = ssym("x")
    A[:,1] = x
    A[1,:] = 5


    print "A = ", NP.dot(A.T,A)
    
  def test_vertcat(self):
    self.message("vertcat")
    A = DMatrix(2,3,1)
    B = DMatrix(4,3)
    C = vertcat([A,B])
    
    self.checkarray(C.shape,(6,3),"vertcat shape")
    self.assertEqual(C.size(),A.size(),"vertcat size")
    
    self.assertRaises(RuntimeError,lambda : horzcat([A,B]))
    
  def test_horzcat(self):
    self.message("horcat")
    A = DMatrix(3,2,1)
    B = DMatrix(3,4)
    C = horzcat([A,B])
    
    self.checkarray(C.shape,(3,6),"horzcat shape")
    self.assertEqual(C.size(),A.size(),"vertcat size")
    
    self.assertRaises(RuntimeError,lambda : vertcat([A,B]))
    
    
  def test_veccat(self):
    self.message("vecccat")
    A = DMatrix(2,3)
    A[0,1] = 2
    A[1,0] = 1
    A[1,2] = 3
    B = DMatrix(3,1)
    B[0,0] = 4
    B[1,0] = 5
    B[2,0] = 6
    C = veccat([A,B])
    
    self.checkarray(C.shape,(9,1),"veccat shape")
    self.assertEqual(C.size(),A.size()+B.size(),"veccat size")
    
    self.checkarray(tuple(C.data()),tuple(arange(1,7)),"numbers shape")
    
  def test_vecNZcat(self):
    self.message("vecNZcat")
    A = DMatrix(2,3)
    A[0,1] = 2
    A[1,0] = 1
    A[1,2] = 3
    B = DMatrix(3,1)
    B[0,0] = 4
    B[1,0] = 5
    B[2,0] = 6
    C = vecNZcat([A,B])
    
    self.checkarray(C.shape,(6,1),"vecNZcat shape")
    self.assertEqual(C.size(),A.size()+B.size(),"vecNZcat size")
    
    self.checkarray(tuple(C.data()),tuple(arange(1,7)),"numbers shape")
    
  def test_IMatrix_indexing(self):
    self.message("IMatrix")
    A = IMatrix(2,2)
    A[0,0] = 1
    A[1,1] = 3
    A[0,1] = 2
    A[1,0] = 4
    
    
    B = DMatrix([1,2,3,4,5])
    
    B_ = B[A]
    
    self.checkarray(B_,DMatrix([[2,3],[5,4]]),"Imatrix indexing")

    B[A] = DMatrix([[1,2],[3,4]])
    
    self.checkarray(B,DMatrix([1,1,2,4,3]),"Imatrix indexing assignement")
    
    #B[A].set(DMatrix([[10,20],[30,40]]))
    
    #self.checkarray(B,DMatrix([1,10,20,40,30]),"Imatrix indexing setting")
    
    B = IMatrix([1,2,3,4,5])
    
    B_ = B[A]
    
    self.checkarray(array(B_),DMatrix([[2,3],[5,4]]),"Imatrix indexing")

    B[A] = IMatrix([[1,2],[3,4]])
    
    self.checkarray(array(B),DMatrix([1,1,2,4,3]),"Imatrix indexing assignement")
  
  def test_IMatrix_index_slice(self):
    self.message("IMatrix combined with slice")

    A = IMatrix(2,2)
    A[0,0] = 0
    A[1,1] = 1
    A[0,1] = 2
    A[1,0] = 0
    
    
    B = DMatrix([[1,2,3],[4,5,6],[7,8,9],[10,11,12]])
    F = DMatrix([[1,2],[4,5]])

    self.checkarray(B[:,A],DMatrix([[1,3],[1,2],[4,6],[4,5],[7,9],[7,8],[10,12],[10,11]]),"B[:,A]")
    self.checkarray(B[A,:],DMatrix([[1,7,2,8,3,9],[1,4,2,5,3,6]]),"B[A,:]")
    

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
    
    C = DMatrix([[1,2,3],[4,5,6],[7,8,9],[10,11,12]])
    print C
    print A.data()
    print B.data()
    print C[list(A.data()),list(B.data())]
    print C[A,B]
    self.checkarray(C[A,B],DMatrix([[3,7],[0,5]]),"C[A,B]")

  def test_issue298(self):
    self.message("Issue #298")
    a = DMatrix(4,1)
    b = c.reshape(a,2,2)
    self.assertEqual(type(a),type(b))

    a = IMatrix(4,1)
    b = DMatrix(a)
    self.assertTrue(isinstance(b,DMatrix))
    
    a = DMatrix(4,1)
    self.assertRaises(RuntimeError,lambda : IMatrix(a))
    
if __name__ == '__main__':
    unittest.main()

