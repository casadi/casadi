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
    
  def test_trans(self):
    self.message("trans")
    a = DMatrix(0,1)
    b = trans(a)
    self.assertEquals(b.size1(),1)
    self.assertEquals(b.size2(),0)
    
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

  def test_nonmonotonous_indexing(self):
    self.message("non-monotonous indexing")
    # Regression test for #354
    A = DMatrix([[1,2,3],[4,5,6],[7,8,9]])
    B = A[[0,2,1],0]
    self.checkarray(DMatrix([1,7,4]),B,"non-monotonous")
    
    B = A[0,[0,2,1]]
    self.checkarray(DMatrix([1,3,2]).T,B,"non-monotonous")
    
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
    
    B = DMatrix(5,1)
   
    self.assertRaises(Exception, lambda : B[A])
  
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
    
    self.assertRaises(Exception, lambda : F[:,A])
    
    self.checkarray(B[A,1],DMatrix([[2,8],[2,5]]),"B[A,1]")
    
    self.checkarray(B[1,A],DMatrix([[4,6],[4,5]]),"B[1,A]")

  def test_IMatrix_index_slice_assignment(self):
    self.message("IMatrix combined with slice assignment")

    A = IMatrix(2,2)
    A[0,0] = 0
    A[1,1] = 1
    A[0,1] = 2
    A[1,0] = 0
    
    
    B = DMatrix([[1,2,3],[4,5,6],[7,8,9],[10,11,12]])
    B_ = DMatrix(B)
    B[:,A] = DMatrix([[1,3],[1,2],[4,6],[4,5],[7,9],[7,8],[10,12],[10,11]])*2
    
    self.checkarray(B,2*B_,"B[:,A] = ")
    
    B[:,A] = 7
    
    self.checkarray(B,DMatrix([[7,7,7],[7,7,7],[7,7,7],[7,7,7]]),"B[:,A] = ")
    
    B[A,:] = DMatrix([[1,7,2,8,3,9],[1,4,2,5,3,6]])
    
    self.checkarray(B,DMatrix([[1,2,3],[4,5,6],[7,8,9],[7,7,7]]),"B[A,:] = ")

    B[A,:] = 6
    
    self.checkarray(B,DMatrix([[6,6,6],[6,6,6],[6,6,6],[7,7,7]]),"B[A,:] = ")
    
    B=DMatrix(3,4)
    B[:,A] = 7
    
    self.checkarray(B,DMatrix([[7,7,7,0],[7,7,7,0],[7,7,7,0]]),"B[:,A] = ")
    
    B=DMatrix(4,4)
    B[A,:] = 8

    self.checkarray(B,DMatrix([[8,8,8,8],[8,8,8,8],[8,8,8,8],[0,0,0,0]]),"B[A,:] = ")
    
    
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
    F = DMatrix([[1,2],[4,5]])

    self.checkarray(C[A,B],DMatrix([[3,7],[0,5]]),"C[A,B]")
    self.assertRaises(Exception, lambda : F[A,B])
    
    C = DMatrix(3,4)
    C_ = C[A,B]
    self.assertEqual(C_.size(),3)
    self.checkarray(C_,DMatrix([[0,0],[0,0]]),"C[A,B]")

  def test_IMatrix_IMatrix_index_assignment(self):
    self.message("IMatrix IMatrix index assignment")

    A = IMatrix(2,2)
    A[0,0] = 0
    A[1,1] = 1
    A[0,1] = 2
    
    B = IMatrix(2,2)
    B[0,0] = 2
    B[1,1] = 1
    B[0,1] = 0
    
    C = DMatrix.zeros((3,4))
    C_ = DMatrix(2,2)
    C_[0,0] = 3
    C_[0,1] = 7
    C_[1,1] = 5;
    
    C[A,B] = C_

    self.checkarray(C[A,B],DMatrix([[3,7],[0,5]]),"C[A,B]")
    
    C = DMatrix(3,4)
    C[A,B] = C_
    self.checkarray(C[A,B],DMatrix([[3,7],[0,5]]),"C[A,B]")

  def test_index_setting(self):
    self.message("index setting")
    B = DMatrix([1,2,3,4,5])
    
    B[0] = 8
    self.checkarray(B,DMatrix([8,2,3,4,5]),"index setting")
    B[1,0] = 4
    self.checkarray(B,DMatrix([8,4,3,4,5]),"index setting")
    B[:,0] = 7
    self.checkarray(B,DMatrix([7,7,7,7,7]),"index setting")
    #B[0].set(3)
    #self.checkarray(B,DMatrix([3,7,7,7,7]),"index setting")
    #B[0].setAll(4)
    #self.checkarray(B,DMatrix([4,4,4,4,4]),"index setting")

    
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
    
  def test_det(self):
    self.message("Determinant")
    npy_det = numpy.linalg.det
    
    a = DMatrix(1,1)
    a[0,0] = 5
    self.checkarray(det(a),npy_det(a),"det()")

    a = DMatrix(5,5)
    for i in range(5):
      a[i,i] = i+1

    self.checkarray(det(a),npy_det(a),"det()")
    
    a = DMatrix(5,5)
    for i in range(4):
      a[i,i] = i+1
    a[0,4] = 3
    a[4,0] = 7
    
    self.checkarray(det(a),npy_det(a),"det()")
    
    a = DMatrix(5,5)
    for i in range(5):
      for j in range(5):
        a[i,j] = i+j
    
    self.checkarray(det(a),npy_det(a),"det()")

    a = DMatrix(5,5)
    for i in range(4):
      for j in range(5):
        a[i,j] = i+j
    
    self.checkarray(det(a),npy_det(a),"det()")
    
    a = DMatrix(5,5)
    for i in range(5):
      for j in range(4):
        a[i,j] = i+j
    
    self.checkarray(det(a),npy_det(a),"det()")
    
    a = DMatrix(5,5)
    for i in range(4):
      for j in range(5):
        a[i,j] = i+j
    a[4,1] = 12
    
    self.checkarray(det(a),npy_det(a),"det()")
    
    a = DMatrix(5,5)
    for i in range(5):
      for j in range(4):
        a[i,j] = i+j
    a[1,4] = 12
    
    self.checkarray(det(a),npy_det(a),"det()")
    
    a = DMatrix(5,5)
    for i in range(4):
      for j in range(5):
        a[i,j] = i+j
    a[4,2] = 12
    
    self.checkarray(det(a),npy_det(a),"det()")
    
    a = DMatrix(5,5)
    for i in range(5):
      for j in range(4):
        a[i,j] = i+j
    a[2,4] = 12
    
    self.checkarray(det(a),npy_det(a),"det()")
    
    a = DMatrix(50,50)
    for i in range(50):
      a[i,i] = i+1

    self.checkarray(det(a)/npy_det(a),1,"det()")
    
  def test_inv_sparsity(self):
    self.message("sparsity pattern of inverse")

    n = 8

    sp = sp_tril(n)

    x  = SXMatrix(sp,[SX("a%d" % i) for i in range(sp.size())])

    
    x_ = DMatrix(x.sparsity(),1)
    
    I_ = DMatrix(inv(x).sparsity(),1)
    
    # For a reducible matrix, struct(A^(-1)) = struct(A) 
    self.checkarray(x_,I_,"inv")
    
    sp = sp_tril(n)

    x  = SXMatrix(sp,[SX("a%d" % i) for i in range(sp.size())])
    x[0,n-1] = 1 
    
    
    I_ = DMatrix(inv(x).sparsity(),1)
    
    # An irreducible matrix has a dense inverse in general
    self.checkarray(DMatrix.ones(n,n),I_,"inv")

    x  = SXMatrix(sp,[SX("a%d" % i) for i in range(sp.size())])
    x[0,n/2] = 1 
    
    s_ = DMatrix(sp,1)
    s_[:,:n/2+1] = 1
    
    I_ = DMatrix(inv(x).sparsity(),1)
    
    makeDense(s_)
    makeDense(I_)
    # An irreducible matrix does not have to be dense per se
    self.checkarray(s_,I_,"inv")

  def test_Imatrix_operations(self):
    self.message("IMatrix operations")
    a = IMatrix(2,2,1)
    b = horzcat([a,a])
    self.assertTrue(isinstance(b,IMatrix))
    
  def test_mul(self):
    A = DMatrix.ones((4,3))
    B = DMatrix.ones((3,8))
    C = DMatrix.ones((8,7))
    
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
    
  def test_remove(self):
    self.message("remove")
    B = DMatrix([[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16],[17,18,19,20]])

    A = DMatrix(B)
    A.remove([],[])
    self.checkarray(A, B,"remove nothing")
    
    A = DMatrix(B)
    A.remove([],[1])
    self.checkarray(A, DMatrix([[1,3,4],[5,7,8],[9,11,12],[13,15,16],[17,19,20]]),"remove a column")
   
    A = DMatrix(B)
    A.remove([0,3],[1])
    self.checkarray(A, DMatrix([[5,7,8],[9,11,12],[17,19,20]]),"remove a column and two rows ")
    
  def test_comparisons(self):
    for m in [DMatrix,IMatrix]:
      A = m([[5,4],[2,1]])
      
      for c in [6,6.0,DMatrix([6]),IMatrix([6]),matrix(6)]:
        self.checkarray(A<=c,m([[1,1],[1,1]]),"<=")
        self.checkarray(A<c,m([[1,1],[1,1]]),"<")
        self.checkarray(A>c,m([[0,0],[0,0]]),">")
        self.checkarray(A>=c,m([[0,0],[0,0]]),">=")
        self.checkarray(A==c,m([[0,0],[0,0]]),"==")
        self.checkarray(A!=c,m([[1,1],[1,1]]),"!=")
        
        self.checkarray(c>=A,m([[1,1],[1,1]]),"<=")
        self.checkarray(c>A,m([[1,1],[1,1]]),"<")
        self.checkarray(c<A,m([[0,0],[0,0]]),">")
        self.checkarray(c<=A,m([[0,0],[0,0]]),">=")
        self.checkarray(c==A,m([[0,0],[0,0]]),"==")
        self.checkarray(c!=A,m([[1,1],[1,1]]),"!=")
        
      for c in [5,5.0,DMatrix([5]),IMatrix([5]),matrix(5)]:
        self.checkarray(A<=c,m([[1,1],[1,1]]),"<=")
        self.checkarray(A<c,m([[0,1],[1,1]]),"<")
        self.checkarray(A>c,m([[0,0],[0,0]]),">")
        self.checkarray(A>=c,m([[1,0],[0,0]]),">=")
        self.checkarray(A==c,m([[1,0],[0,0]]),"==")
        self.checkarray(A!=c,m([[0,1],[1,1]]),"!=")

        self.checkarray(c>=A,m([[1,1],[1,1]]),"<=")
        self.checkarray(c>A,m([[0,1],[1,1]]),"<")
        self.checkarray(c<A,m([[0,0],[0,0]]),">")
        self.checkarray(c<=A,m([[1,0],[0,0]]),">=")
        self.checkarray(c==A,m([[1,0],[0,0]]),"==")
        self.checkarray(c!=A,m([[0,1],[1,1]]),"!=")
        
      for c in [4,4.0,DMatrix([4]),IMatrix([4]),matrix(4)]:
        self.checkarray(A<=c,m([[0,1],[1,1]]),"<=")
        self.checkarray(A<c,m([[0,0],[1,1]]),"<")
        self.checkarray(A>c,m([[1,0],[0,0]]),">")
        self.checkarray(A>=c,m([[1,1],[0,0]]),">=")
        self.checkarray(A==c,m([[0,1],[0,0]]),"==")
        self.checkarray(A!=c,m([[1,0],[1,1]]),"!=")

        self.checkarray(c>=A,m([[0,1],[1,1]]),"<=")
        self.checkarray(c>A,m([[0,0],[1,1]]),"<")
        self.checkarray(c<A,m([[1,0],[0,0]]),">")
        self.checkarray(c<=A,m([[1,1],[0,0]]),">=")
        self.checkarray(c==A,m([[0,1],[0,0]]),"==")
        self.checkarray(c!=A,m([[1,0],[1,1]]),"!=")
        
      for c in [1,1.0,DMatrix([1]),IMatrix([1]),matrix(1)]:
        self.checkarray(A<=c,m([[0,0],[0,1]]),"<=")
        self.checkarray(A<c,m([[0,0],[0,0]]),"<")
        self.checkarray(A>c,m([[1,1],[1,0]]),">")
        self.checkarray(A>=c,m([[1,1],[1,1]]),">=")
        self.checkarray(A==c,m([[0,0],[0,1]]),"==")
        self.checkarray(A!=c,m([[1,1],[1,0]]),"!=")

        self.checkarray(c>=A,m([[0,0],[0,1]]),"<=")
        self.checkarray(c>A,m([[0,0],[0,0]]),"<")
        self.checkarray(c<A,m([[1,1],[1,0]]),">")
        self.checkarray(c<=A,m([[1,1],[1,1]]),">=")
        self.checkarray(c==A,m([[0,0],[0,1]]),"==")
        self.checkarray(c!=A,m([[1,1],[1,0]]),"!=")
        
      for c in [0,DMatrix([0]),IMatrix([0]),matrix(0)]:
        self.checkarray(A<=c,m([[0,0],[0,0]]),"<=")
        self.checkarray(A<c,m([[0,0],[0,0]]),"<")
        self.checkarray(A>c,m([[1,1],[1,1]]),">")
        self.checkarray(A>=c,m([[1,1],[1,1]]),">=")
        self.checkarray(A==c,m([[0,0],[0,0]]),"==")
        self.checkarray(A!=c,m([[1,1],[1,1]]),"!=")

        self.checkarray(c>=A,m([[0,0],[0,0]]),"<=")
        self.checkarray(c>A,m([[0,0],[0,0]]),"<")
        self.checkarray(c<A,m([[1,1],[1,1]]),">")
        self.checkarray(c<=A,m([[1,1],[1,1]]),">=")
        self.checkarray(c==A,m([[0,0],[0,0]]),"==")
        self.checkarray(c!=A,m([[1,1],[1,1]]),"!=")

  def test_all_any(self):
    for m in [DMatrix,IMatrix]:
      A = m([[1,1],[1,1]])
      self.assertTrue(all(A))
      self.assertTrue(any(A))
      
      A = m([[1,0],[1,1]])
      self.assertFalse(all(A))
      self.assertTrue(any(A))

      A = m([[0,0],[0,0]])
      self.assertFalse(all(A))
      self.assertFalse(any(A))
    
if __name__ == '__main__':
    unittest.main()

