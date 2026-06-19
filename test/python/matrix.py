#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             KU Leuven. All rights reserved.
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
from casadi import mtimes
import casadi as ca
import numpy as np
from numpy import inf, pi
import casadi as c
import numpy
from numpy import eye, linalg, arange, matrix
import unittest
from types import *
from helpers import *
import numpy
from itertools import *

class Matrixtests(casadiTestCase):
  def test_constructorlol(self):
    self.message("List of list constructor")
    a=ca.DM(array([[1,2,3],[4,5,6],[7,8,9]]))
    b=ca.DM([[1,2,3],[4,5,6],[7,8,9]])
    self.checkarray(a,b,"List of list constructor")

  def test_sum(self):
    self.message("sum")
    D=ca.DM([[1,2,3],[4,5,6],[7,8,9]])
    self.checkarray(c.sum1(D),array([[12,15,18]]),'sum()')
    self.checkarray(c.sum2(D),array([[6,15,24]]).T,'sum()')


  def test_inv(self):
    self.message("Matrix inverse")
    a = ca.DM([[1,2],[1,3]])
    self.checkarray(c.inv(a) @ a,eye(2),"DM inverse")

  def test_trans(self):
    self.message("trans")
    a = ca.DM(0,1)
    b = a.T
    self.assertEqual(b.size1(),1)
    self.assertEqual(b.size2(),0)

  def test_vertcat(self):
    self.message("vertcat")
    A = ca.DM.ones(2,3)
    B = ca.DM(4,3)
    C = ca.vertcat(A,B)

    self.checkarray(C.shape,(6,3),"vertcat shape")
    self.assertEqual(C.nnz(),A.nnz(),"vertcat size")

    self.assertRaises(RuntimeError,lambda : ca.horzcat(A,B))

  def test_horzcat(self):
    self.message("horcat")
    A = ca.DM.ones(3,2)
    B = ca.DM(3,4)
    C = ca.horzcat(A,B)

    self.checkarray(C.shape,(3,6),"horzcat shape")
    self.assertEqual(C.nnz(),A.nnz(),"vertcat size")

    self.assertRaises(RuntimeError,lambda : ca.vertcat(A,B))

  def test_diagcat(self):

    x = ca.MX.sym("x",2,2)
    y = ca.MX.sym("y",ca.Sparsity.lower(3))
    z = ca.MX.sym("z",4,2)

    L = [x,y,z]

    fMX = ca.Function("fMX", L,[ca.diagcat(*L)])

    LSX = [ ca.SX.sym("",i.sparsity()) for i in L ]
    fSX = ca.Function("fSX", LSX,[ca.diagcat(*LSX)])

    for f in [fMX,fSX]:
      for i in range(3):
        f_in[i]=list(range(f.nnz_in(i)))  # pyright: ignore[reportUndefinedVariable]

    self.checkfunction(fMX,fSX)

  def test_veccat(self):
    self.message("vecccat")
    A = ca.DM(2,3)
    A[0,1] = 2
    A[1,0] = 1
    A[1,2] = 3
    B = ca.DM(3,1)
    B[0,0] = 4
    B[1,0] = 5
    B[2,0] = 6
    C = ca.veccat(*[A,B])

    self.checkarray(C.shape,(9,1),"veccat shape")
    self.assertEqual(C.nnz(),A.nnz()+B.nnz(),"veccat size")

    self.checkarray(tuple(C.nonzeros()),tuple(arange(1,7)),"numbers shape")

  def test_slicestepnegative(self):
    self.message("Slice step negative")
    a1 = [1,2,3,4,5]
    a2 = ca.DM(a1)

    self.checkarray(a2[0:4:-1,0],ca.DM(a1[0:4:-1])) # gives empty set
    self.checkarray(a2[4:0:-1,0],ca.DM(a1[4:0:-1])) # gives [5, 4, 3, 2]
    self.checkarray(a2[0:4:-2,0],ca.DM(a1[0:4:-2])) # gives empty set
    self.checkarray(a2[4:0:-2,0],ca.DM(a1[4:0:-2])) # gives [5, 4, 3, 2]
    self.checkarray(a2[1:4:-2,0],ca.DM(a1[1:4:-2])) # gives empty set
    self.checkarray(a2[4:1:-2,0],ca.DM(a1[4:1:-2])) # gives [5, 4, 3, 2]
    self.checkarray(a2[0:3:-2,0],ca.DM(a1[0:3:-2])) # gives empty set
    self.checkarray(a2[3:0:-2,0],ca.DM(a1[3:0:-2])) # gives [5, 4, 3, 2]
    self.checkarray(a2[::-1,0],ca.DM(a1[::-1])) # gives [5, 4, 3, 2, 1]
    self.checkarray(a2[::1,0],ca.DM(a1[::1])) # gives [1,2,3,4,5]
    self.checkarray(a2[2::-1,0],ca.DM(a1[2::-1])) # gives [3,2,1]
    self.checkarray(a2[:2:-1,0],ca.DM(a1[:2:-1])) # gives [5,4]
    self.checkarray(a2[-1::-1,0],ca.DM(a1[-1::-1])) # gives [5,4,3,2,1]
    self.checkarray(a2[:-1:-1,0],ca.DM(a1[:-1:-1])) # gives []

  def test_indexingOutOfBounds(self):
    self.message("Indexing out of bounds")
    y = ca.DM.zeros(4, 5)
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

  def huge_slice(self):
    self.message("huge slice")
    a = ca.SX.sym("a",ca.Sparsity.diag(50000))

    a[:,:]

  def test_nonmonotonous_indexing(self):
    self.message("non-monotonous indexing")
    # Regression test for #354
    A = ca.DM([[1,2,3],[4,5,6],[7,8,9]])
    B = A[[0,2,1],0]
    self.checkarray(ca.DM([1,7,4]),B,"non-monotonous")

    B = A[0,[0,2,1]]
    self.checkarray(ca.DM([1,3,2]).T,B,"non-monotonous")

  def test_IM_indexing(self):
    self.message("IM")
    A = ca.DM(2,2)
    A[0,0] = 1
    A[1,1] = 3
    A[0,1] = 2
    A[1,0] = 4


    B = ca.DM([1,2,3,4,5])

    B_ = B[A]

    self.checkarray(B_,ca.DM([[2,3],[5,4]]),"matrix indexing")

    B[A] = ca.DM([[1,2],[3,4]])

    self.checkarray(B,ca.DM([1,1,2,4,3]),"matrix indexing assignement")

    #B[A].set(DM([[10,20],[30,40]]))

    #self.checkarray(B,DM([1,10,20,40,30]),"Imatrix indexing setting")

    B = ca.DM([1,2,3,4,5])

    B_ = B[A]

    self.checkarray(array(B_),ca.DM([[2,3],[5,4]]),"Imatrix indexing")

    B[A] = ca.DM([[1,2],[3,4]])

    self.checkarray(array(B),ca.DM([1,1,2,4,3]),"Imatrix indexing assignement")

    B = ca.DM(5,1)

    self.assertRaises(Exception, lambda : B.nz[A])

  def test_sparsity_indexing(self):
    self.message("sparsity")

    B = ca.DM([[1,2,3,4,5],[6,7,8,9,10]])

    A = ca.DM([[1,1,0,0,0],[0,0,1,0,0]])
    A = ca.sparsify(A)
    sp = A.sparsity()


    self.checkarray(B[sp],ca.DM([[1,2,0,0,0],[0,0,8,0,0]]),"sparsity indexing")

    B[sp] = -4

    self.checkarray(B,ca.DM([[-4,-4,3,4,5],[6,7,-4,9,10]]),"sparsity indexing assignement")

    B = ca.DM([[1,2,3,4,5],[6,7,8,9,10]])

    B[sp] = 2*B

    self.checkarray(B,ca.DM([[2,4,3,4,5],[6,7,16,9,10]]),"Imatrix indexing assignement")

    self.assertRaises(Exception, lambda : B[ca.Sparsity.dense(4,4)])

  def test_index_setting(self):
    self.message("index setting")
    B = ca.DM([1,2,3,4,5])

    B[0] = 8
    self.checkarray(B,ca.DM([8,2,3,4,5]),"index setting")
    B[1,0] = 4
    self.checkarray(B,ca.DM([8,4,3,4,5]),"index setting")
    B[:,0] = 7
    self.checkarray(B,ca.DM([7,7,7,7,7]),"index setting")
    #B[0].set(3)
    #self.checkarray(B,DM([3,7,7,7,7]),"index setting")
    #B[0].setAll(4)
    #self.checkarray(B,DM([4,4,4,4,4]),"index setting")


  def test_issue298(self):
    self.message("Issue #298")
    a = ca.DM(4,1)
    b = c.reshape(a,2,2)
    self.assertEqual(type(a),type(b))

  def test_det(self):
    self.message("Determinant")
    npy_det = numpy.linalg.det

    a = ca.DM(1,1)
    a[0,0] = 5
    self.checkarray(ca.det(a),npy_det(a),"det()")

    a = ca.DM(5,5)
    for i in range(5):
      a[i,i] = i+1

    self.checkarray(ca.det(a),npy_det(a),"det()")

    a = ca.DM(5,5)
    for i in range(4):
      a[i,i] = i+1
    a[0,4] = 3
    a[4,0] = 7

    self.checkarray(ca.det(a),npy_det(a),"det()")

    a = ca.DM(5,5)
    for i in range(5):
      for j in range(5):
        a[i,j] = i+j

    self.checkarray(ca.det(a),npy_det(a),"det()")

    a = ca.DM(5,5)
    for i in range(4):
      for j in range(5):
        a[i,j] = i+j

    self.checkarray(ca.det(a),npy_det(a),"det()")

    a = ca.DM(5,5)
    for i in range(5):
      for j in range(4):
        a[i,j] = i+j

    self.checkarray(ca.det(a),npy_det(a),"det()")

    a = ca.DM(5,5)
    for i in range(4):
      for j in range(5):
        a[i,j] = i+j
    a[4,1] = 12

    self.checkarray(ca.det(a),npy_det(a),"det()")

    a = ca.DM(5,5)
    for i in range(5):
      for j in range(4):
        a[i,j] = i+j
    a[1,4] = 12

    self.checkarray(ca.det(a),npy_det(a),"det()")

    a = ca.DM(5,5)
    for i in range(4):
      for j in range(5):
        a[i,j] = i+j
    a[4,2] = 12

    self.checkarray(ca.det(a),npy_det(a),"det()")

    a = ca.DM(5,5)
    for i in range(5):
      for j in range(4):
        a[i,j] = i+j
    a[2,4] = 12

    self.checkarray(ca.det(a),npy_det(a),"det()")

    a = ca.DM(50,50)
    for i in range(50):
      a[i,i] = i+1

    self.checkarray(ca.det(a)/npy_det(a),1,"det()")

  @skip(not ca.GlobalOptions.getSimplificationOnTheFly())
  def test_inv_sparsity(self):
    self.message("sparsity pattern of inverse")

    n = 8

    sp = ca.Sparsity.lower(n)

    x  = ca.SX(sp,ca.vertcat(*[ca.SX.sym("a%d" % i) for i in range(sp.nnz())]))


    x_ = ca.DM.ones(x.sparsity())

    I_ = ca.DM.ones(ca.inv(x).sparsity())

    # For a reducible matrix, struct(A^(-1)) = struct(A)
    self.checkarray(x_,I_,"inv")

    sp = ca.Sparsity.lower(n)

    x = ca.SX.sym("a", sp)
    x[0,n-1] = 1


    I_ = ca.DM.ones(ca.inv(x).sparsity())

    # An irreducible matrix has a dense inverse in general
    self.checkarray(ca.DM.ones(n,n),I_,"inv")

    x = ca.SX.sym("a", sp)
    x[0,int(n/2)] = 1

    s_ = ca.DM.ones(sp)
    s_[:,:int(n/2)+1] = 1

    I_ = ca.DM.ones(ca.inv_minor(x).sparsity())

    s_ = ca.densify(s_)
    T_ = ca.densify(I_)
    # An irreducible matrix does not have to be dense per se
    self.checkarray(s_,I_,"inv")

  def test_mtimes(self):
    A = ca.DM.ones((4,3))
    B = ca.DM.ones((3,8))
    C = ca.DM.ones((8,7))

    self.assertRaises(RuntimeError,lambda : mtimes([]))

    D = mtimes([A])

    self.assertEqual(D.shape[0],4)
    self.assertEqual(D.shape[1],3)

    D = mtimes([A,B])

    self.assertEqual(D.shape[0],4)
    self.assertEqual(D.shape[1],8)

    D = mtimes([A,B,C])

    self.assertEqual(D.shape[0],4)
    self.assertEqual(D.shape[1],7)

  def test_remove(self):
    self.message("remove")
    B = ca.DM([[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16],[17,18,19,20]])

    A = ca.DM(B)
    A.remove([],[])
    self.checkarray(A, B,"remove nothing")

    A = ca.DM(B)
    A.remove([],[1])
    self.checkarray(A, ca.DM([[1,3,4],[5,7,8],[9,11,12],[13,15,16],[17,19,20]]),"remove a column")

    A = ca.DM(B)
    A.remove([0,3],[1])
    self.checkarray(A, ca.DM([[5,7,8],[9,11,12],[17,19,20]]),"remove a column and two rows ")

  def test_comparisons(self):
    m = ca.DM
    A = m([[5,4],[2,1]])

    for c in [6,6.0,ca.DM([6]),np.array([6]),matrix(6)]:
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

    for c in [5,5.0,ca.DM([5]),np.array([5])]+([matrix(5)] if check_matrix else []):
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

    for c in [4,4.0,ca.DM([4]),np.array([4]),matrix(4)]:
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

    for c in [1,1.0,ca.DM([1]),np.array([1]),matrix(1)]:
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

    for c in [0,ca.DM([0]),np.array([0]),matrix(0)]:
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

  def test_truth(self):
    self.assertTrue(bool(ca.DM([1])))
    self.assertFalse(bool(ca.DM([0])))
    self.assertTrue(bool(ca.DM([0.2])))
    self.assertTrue(bool(ca.DM([-0.2])))
    self.assertRaises(Exception, lambda : bool(ca.DM([2.0,3])))
    self.assertRaises(Exception, lambda : bool(ca.DM()))

  def test_listslice(self):
    def check(d,rowbase,colbase):
      for col in permutations(colbase):
        for row in permutations(rowbase):
          r = ca.DM.zeros(len(row),len(col))
          for i,ii in enumerate(row):
            for j,jj in enumerate(col):
              r[i,j] = d[ii,jj]
          self.checkarray(d[row,col],r,"%s[%s,%s]" % (repr(d),str(row),str(col)))


    # get1
    check(ca.DM(ca.Sparsity.dense(3,3),list(range(3*3))),[0,1,2],[0,1,2])
    check(ca.DM(ca.Sparsity.dense(4,4),list(range(4*4))),[0,1,3],[0,2,3])
    check(ca.DM(ca.Sparsity.dense(3,3),list(range(3*3))),[0,0,1],[0,0,1])
    check(ca.DM(ca.Sparsity.dense(3,3),list(range(3*3))),[0,0,2],[0,0,2])
    check(ca.DM(ca.Sparsity.dense(3,3),list(range(3*3))),[1,1,2],[1,1,2])

    sp = ca.Sparsity.lower(4)
    d = ca.DM(sp,list(range(sp.nnz())))
    check(d,[0,1,3],[0,2,3])
    check(d.T,[0,1,3],[0,2,3])

    sp = ca.Sparsity.rowcol([0,1,2],[0,1],4,4)
    d = ca.DM(sp,list(range(sp.nnz())))
    check(d,[0,3],[0,2])

    # get2
    check(ca.DM(ca.Sparsity.dense(2,2),list(range(2*2))),[0,0,0],[0,0,0])
    check(ca.DM(ca.Sparsity.dense(2,2),list(range(2*2))),[0,0,1],[0,0,1])
    check(ca.DM(ca.Sparsity.dense(2,2),list(range(2*2))),[1,1,0],[1,1,0])
    check(ca.DM(ca.Sparsity.dense(2,2),list(range(2*2))),[1,1,1],[1,1,1])

    sp = ca.Sparsity.lower(3)
    d = ca.DM(sp,list(range(sp.nnz())))
    check(d,[0,1,2],[0,1,2])
    check(d.T,[0,1,2],[0,1,2])

    sp = ca.Sparsity.rowcol([0,2],[0,1],4,4)
    d = ca.DM(sp,list(range(sp.nnz())))
    check(d,[0,1,3],[0,2,3])

  def test_sparsesym(self):
    # feature removed in 73f407e
    return
    self.message("sparsesym")
    D = ca.DM([[1,2,-3],[2,-1,0],[-3,0,5]])
    D = ca.sparsify(D)

    i = D.getSym()
    #self.checkarray(list(i),[1,2,-1,-3,5])
    A = 2*D
    A.setSym(i)
    self.checkarray(A,D)

  def test_diagcat(self):
    self.message("diagcat")
    C = ca.diagcat(*[ca.DM([[-1.4,-3.2],[-3.2,-28]]),ca.DM([[15,-12,2.1],[-12,16,-3.8],[2.1,-3.8,15]]),1.8,-4.0])
    r = ca.DM([[-1.4,-3.2,0,0,0,0,0],[-3.2,-28,0,0,0,0,0],[0,0,15,-12,2.1,0,0],[0,0,-12,16,-3.8,0,0],[0,0,2.1,-3.8,15,0,0],[0,0,0,0,0,1.8,0],[0,0,0,0,0,0,-4]])
    r = ca.sparsify(r)
    self.checkarray(C,r)

  def test_diag_sparse(self):
    self.message("diag sparse")

    for n in [[0,1,0,0,2,3,4,5,6,0],[1,2,3,0],[0,1,2,3]]:
      d = ca.DM(n)
      D = ca.DM(n)
      d = ca.sparsify(d)
      m = c.diag(d)
      M = ca.sparsify(c.diag(D))

      self.checkarray(m.sparsity().colind(),M.sparsity().colind())
      self.checkarray(m.sparsity().row(),M.sparsity().row())

  def test_sprank(self):
    self.message("sprank")

    a = ca.DM([[1,0,0],[0,1,0],[0,0,1]])
    a = ca.sparsify(a)
    self.assertEqual(ca.sprank(a),3)

    a = ca.DM([[1,0,0],[0,0,0],[0,0,1]])
    a = ca.sparsify(a)
    self.assertEqual(ca.sprank(a),2)

    a = ca.DM([[0,0,0],[0,0,0],[0,0,1]])
    a = ca.sparsify(a)
    self.assertEqual(ca.sprank(a),1)

    a = ca.DM([[0,0,0],[0,0,0],[0,0,0]])
    a = ca.sparsify(a)
    self.assertEqual(ca.sprank(a),0)

    self.assertEqual(ca.sprank(ca.DM.ones(1,3)),1)
    self.assertEqual(ca.sprank(ca.DM.ones(3,1)),1)
    self.assertEqual(ca.sprank(ca.DM.ones(2,3)),2)
    self.assertEqual(ca.sprank(ca.DM.ones(3,2)),2)
    self.assertEqual(ca.sprank(ca.DM.ones(3,3)),3)
    self.assertEqual(ca.sprank(ca.DM.ones(3,3)),3)

    A = ca.DM(6,4)
    A[0,0] = 1
    A[1,2] = 1
    A[2,2] = 1
    A[5,3] = 1

    self.assertEqual(ca.sprank(A),3)

  def test_cross(self):
    self.message("cross products")

    crossc = c.cross

    self.checkarray(crossc(ca.DM([1,0,0]),ca.DM([0,1,0])),ca.DM([0,0,1]))

    self.checkarray(crossc(ca.DM([1.1,1.3,1.7]),ca.DM([2,3,13])),ca.DM([11.8,-10.9,0.7]))
    self.checkarray(crossc(ca.DM([1.1,1.3,1.7]).T,ca.DM([2,3,13]).T),ca.DM([11.8,-10.9,0.7]).T)

    self.checkarray(crossc(ca.DM([[1.1,1.3,1.7],[1,0,0],[0,0,1],[4,5,6]]),ca.DM([[2,3,13],[0,1,0],[0,0,1],[1,0,1]])),ca.DM([[11.8,-10.9,0.7],[0,0,1],[0,0,0],[5,2,-5]]))
    self.checkarray(crossc(ca.DM([[1.1,1.3,1.7],[1,0,0],[0,0,1],[4,5,6]]).T,ca.DM([[2,3,13],[0,1,0],[0,0,1],[1,0,1]]).T),ca.DM([[11.8,-10.9,0.7],[0,0,1],[0,0,0],[5,2,-5]]).T)

    self.checkarray(crossc(ca.DM([[1.1,1.3,1.7],[1,0,0],[0,0,1],[4,5,6]]),ca.DM([[2,3,13],[0,1,0],[0,0,1],[1,0,1]]),2),ca.DM([[11.8,-10.9,0.7],[0,0,1],[0,0,0],[5,2,-5]]))

    self.checkarray(crossc(ca.DM([[1.1,1.3,1.7],[1,0,0],[0,0,1],[4,5,6]]).T,ca.DM([[2,3,13],[0,1,0],[0,0,1],[1,0,1]]).T,1),ca.DM([[11.8,-10.9,0.7],[0,0,1],[0,0,0],[5,2,-5]]).T)

  def test_is_regular(self):
    self.assertTrue(ca.DM([1,2]).is_regular())
    self.assertFalse(ca.DM([1,inf]).is_regular())
    self.assertFalse(ca.DM.nan(2).is_regular())

  def test_sizes(self):
    self.assertEqual(ca.Sparsity.diag(10).nnz_diag(),10)
    self.assertEqual(ca.Sparsity.diag(10).nnz_upper(),10)
    self.assertEqual(ca.Sparsity.diag(10).nnz_lower(),10)
    self.assertEqual(ca.Sparsity.dense(10,10).nnz_lower(),10*11/2)
    self.assertEqual(ca.Sparsity.dense(10,10).nnz_upper(),10*11/2)
    self.assertEqual(ca.Sparsity.dense(10,10).nnz_diag(),10)

    self.assertEqual(ca.sparsify(ca.DM([[1,1,0],[1,0,1],[0,0,0]])).nnz_diag(),1)
    self.assertEqual(ca.sparsify(ca.DM([[1,1,0],[1,0,1],[0,0,0]])).nnz_lower(),2)
    self.assertEqual(ca.sparsify(ca.DM([[1,1,0],[1,0,1],[0,0,0]])).nnz_upper(),3)

  def test_tril2symm(self):
    a = ca.DM(ca.Sparsity.upper(3),list(range(ca.Sparsity.upper(3).nnz()))).T
    s = ca.tril2symm(a)
    self.checkarray(s,ca.DM([[0,1,3],[1,2,4],[3,4,5]]))

    with self.assertRaises(Exception):
      ca.tril2symm(ca.DM.ones(5,3))

    print(ca.DM.ones(5,5).nnz_upper()-ca.DM.ones(5,5).nnz_diag())

    with self.assertRaises(Exception):
      ca.tril2symm(ca.DM.ones(5,5))

  def test_not_null(self):
    x = ca.MX.sym('x',3,1)
    sp = ca.Sparsity.upper(2)
    ca.MX(sp,x)

  def test_segfault(self):
    x = ca.MX.sym('x',10,1)
    sp = ca.Sparsity.upper(2)
    y = ca.triu2symm(ca.MX(sp,x[1:4]))
    f = ca.Function("f", [x],[y])

  def test_append_empty(self):
    a = ca.vertcat(ca.DM(0,0),ca.DM(0,2))

    self.assertEqual(a.size1(),0)
    self.assertEqual(a.size2(),2)

    a = ca.vertcat(ca.DM(0,0),ca.DM(2,0),ca.DM(3,0))

    self.assertEqual(a.size1(),5)
    self.assertEqual(a.size2(),0)

  def test_vertcat_empty(self):
    a = ca.DM(0,2)
    v = ca.vertcat(a,a)

    self.assertEqual(v.size1(),0)
    self.assertEqual(v.size2(),2)

    a = ca.DM(2,0)
    v = ca.vertcat(a,a)

    self.assertEqual(v.size1(),4)
    self.assertEqual(v.size2(),0)

  def test_vertsplit(self):
    a = ca.DM(ca.Sparsity.upper(5),list(range(int(5*6/2)))).T
    v = ca.vertsplit(a,[0,2,4,5])

    self.assertEqual(len(v),3)
    self.checkarray(v[0],ca.DM([[0,0,0,0,0],[1,2,0,0,0]]))
    self.checkarray(v[1],ca.DM([[3,4,5,0,0],[6,7,8,9,0]]))
    self.checkarray(v[2],ca.DM([[10,11,12,13,14]]))

    v = ca.vertsplit(a)
    self.assertEqual(len(v),a.size1())
    self.checkarray(v[0],ca.DM([[0,0,0,0,0]]))
    self.checkarray(v[1],ca.DM([[1,2,0,0,0]]))
    self.checkarray(v[2],ca.DM([[3,4,5,0,0]]))
    self.checkarray(v[3],ca.DM([[6,7,8,9,0]]))
    self.checkarray(v[4],ca.DM([[10,11,12,13,14]]))

    v = ca.vertsplit(a,[0,2,4,5])
    self.assertEqual(len(v),3)
    self.checkarray(v[0],ca.DM([[0,0,0,0,0],[1,2,0,0,0]]))
    self.checkarray(v[1],ca.DM([[3,4,5,0,0],[6,7,8,9,0]]))
    self.checkarray(v[2],ca.DM([[10,11,12,13,14]]))

    v = ca.vertsplit(a,[0,0,3,5])
    self.assertEqual(len(v),3)
    self.assertEqual(v[0].size1(),0)
    self.assertEqual(v[0].size2(),5)
    self.checkarray(v[1],ca.DM([[0,0,0,0,0],[1,2,0,0,0],[3,4,5,0,0]]))
    self.checkarray(v[2],ca.DM([[6,7,8,9,0],[10,11,12,13,14]]))

  def test_horzsplit(self):
    a = ca.DM(ca.Sparsity.upper(5),list(range(int(5*6/2)))).T
    v = ca.horzsplit(a,[0,2,4,5])

    self.assertEqual(len(v),3)
    self.checkarray(v[0],ca.DM([[0,0],[1,2],[3,4],[6,7],[10,11]]))
    self.checkarray(v[1],ca.DM([[0,0],[0,0],[5,0],[8,9],[12,13]]))
    self.checkarray(v[2],ca.DM([[0],[0],[0],[0],[14]]))

    v = ca.horzsplit(a)
    self.assertEqual(len(v),a.size1())
    self.checkarray(v[0],ca.DM([0,1,3,6,10]))
    self.checkarray(v[1],ca.DM([0,2,4,7,11]))
    self.checkarray(v[2],ca.DM([0,0,5,8,12]))
    self.checkarray(v[3],ca.DM([0,0,0,9,13]))
    self.checkarray(v[4],ca.DM([0,0,0,0,14]))

    v = ca.horzsplit(a,[0,2,4,5])
    self.assertEqual(len(v),3)
    self.checkarray(v[0],ca.DM([[0,0],[1,2],[3,4],[6,7],[10,11]]))
    self.checkarray(v[1],ca.DM([[0,0],[0,0],[5,0],[8,9],[12,13]]))
    self.checkarray(v[2],ca.DM([[0],[0],[0],[0],[14]]))

    v = ca.horzsplit(a,[0,0,3,5])
    self.assertEqual(len(v),3)
    self.assertEqual(v[0].size1(),5)
    self.assertEqual(v[0].size2(),0)
    self.checkarray(v[1],ca.DM([[0,0,0],[1,2,0],[3,4,5],[6,7,8],[10,11,12]]))
    self.checkarray(v[2],ca.DM([[0,0],[0,0],[0,0],[9,0],[13,14]]))

  def test_blocksplit(self):
    a = ca.DM(ca.Sparsity.upper(5),list(range(int(5*6/2)))).T
    v = ca.blocksplit(a,[0,2,4,5],[0,1,3,5])

    self.checkarray(v[0][0],ca.DM([0,1]))
    self.checkarray(v[0][1],ca.DM([[0,0],[2,0]]))
    self.checkarray(v[1][0],ca.DM([3,6]))
    self.checkarray(ca.blockcat(v),a)

  def test_solve(self):
    import random

    spA = [
      ca.Sparsity.dense(1,1)
    ]

    for n in range(2,5):
      spA+= [
        ca.Sparsity.diag(n),
        ca.Sparsity.dense(n,n),
        ca.Sparsity.lower(n),
        ca.Sparsity.lower(n).T,
        ca.Sparsity.banded(n,1),
        ca.diagcat(*[ca.Sparsity.diag(n),ca.Sparsity.dense(n,n)]),
        ca.diagcat(*[ca.Sparsity.diag(n),ca.Sparsity.lower(n)]),
        ca.diagcat(*[ca.Sparsity.diag(n),ca.Sparsity.lower(n).T]),
        ca.diagcat(*[ca.Sparsity.lower(n),ca.Sparsity.lower(n).T]),
        ca.Sparsity.diag(n)+ca.Sparsity.rowcol([0],[n-1],n,n),
        ca.Sparsity.diag(n)+ca.Sparsity.rowcol([0,n-1],[n-1,0],n,n),
        ca.Sparsity.diag(n)+ca.Sparsity.triplet(n,n,[0],[n-1]),
        ca.Sparsity.diag(n)+ca.Sparsity.triplet(n,n,[0,n-1],[n-1,0]),
      ]

    for sA in spA:

      random.seed(1)
      a = ca.DM(sA,[random.random() for i in range(sA.nnz())])
      A = ca.SX.sym("a",a.sparsity())
      for sB in [ ca.Sparsity.dense(a.size1(),1), ca.vertcat(ca.Sparsity.dense(1,1),ca.Sparsity(a.size1()-1,1)),ca.Sparsity.lower(a.size1()),ca.Sparsity.lower(a.size1()).T]:

        b = ca.DM(sB,[random.random() for i in range(sB.nnz())])
        B = ca.SX.sym("B",b.sparsity())
        C = ca.solve(A,B)

        f = ca.Function("f", [A,B],[C])


        c = f(a,b)

        c_ref = ca.DM(linalg.solve(a,b))
        c_ref = ca.sparsify(c_ref)

        print(sA.dim(), sB.dim())



        #try:
        print("foo",c,c_ref)
        self.checkarray(c,c_ref)
        #self.assertTrue(min((IM.ones(c_ref.sparsity())-IM.ones(c.sparsity())).nonzeros())==0)
        #except Exception as e:
        #  c.print_dense()
        #  print("sol:")
        #  c.sparsity().spy()
        #  print("ref:")
        #  c_ref.sparsity().spy()
        #  c_ref.print_dense()
        #  a.print_dense()
        #  raise e

  def test_kron(self):
    a = ca.sparsify(ca.DM([[1,0,6],[2,7,0]]))
    b = ca.sparsify(ca.DM([[1,0,0],[2,3,7],[0,0,9],[1,12,13]]))

    c_ = c.kron(a,b)

    self.assertEqual(c_.size1(),a.size1()*b.size1())
    self.assertEqual(c_.size2(),a.size2()*b.size2())
    self.assertEqual(c_.nnz(),a.nnz()*b.nnz())

    self.checkarray(c_,numpy.kron(a,b))

  def test_vec_kron(self):
    A = ca.SX.sym("A",2,3)
    B = ca.SX.sym("B",4,5)
    P = ca.SX.sym("P",A.size2(),B.size1())

    f = ca.Function("f", [ca.vec(P.T),A,B],[ca.vec(mtimes([A,P,B]).T)])

    J = jacobian_old(f, 0, 0)
    J_in = []
    J_in.append(numpy.random.rand(*ca.vec(P.T).shape))
    J_in.append(numpy.random.rand(*A.shape))
    J_in.append(numpy.random.rand(*B.shape))

    res, _  =  J(*J_in)

    ref =  ca.kron(J_in[1],J_in[2].T)

    self.checkarray(res,ref)

  def test_repmat(self):
    a = ca.DM([[1,2],[3,4],[5,6]])
    self.checkarray(ca.repmat(a,2,3),ca.kron(ca.DM.ones(2,3),a))

  def test_triu(self):
    a = ca.DM([[1,2],[3,4]])
    b = ca.triu(a)
    self.checkarray(b, ca.DM([[1,2],[0,4]]) )


  def test_tril(self):
    a = ca.DM([[1,2],[3,4]])
    b = ca.tril(a)
    self.checkarray(b, ca.DM([[1,0],[3,4]]) )

  def test_nz(self):
    a = ca.sparsify(ca.DM([[1,2],[0,0],[3,4]]))
    self.checkarray(a.nz[:], ca.DM([1,3,2,4]) )
    self.checkarray(len(a.nz), 4 )
    self.checkarray(a.nz[:-1], ca.DM([1,3,2]) )
    self.checkarray(a.nz[0], ca.DM([1]) )

  def test_norm_inf_mul(self):
    numpy.random.seed(0)

    A = numpy.random.random((10,2))
    B = numpy.random.random((2,8))

    self.checkarray(ca.norm_inf_mul(A,B),ca.norm_inf(A @ B))
    self.checkarray(ca.DM(ca.norm_0_mul(A,B)),mtimes(A,B).nnz())

    # Sparse
    for i in range(5):
      A[numpy.random.randint(A.shape[0]),numpy.random.randint(A.shape[1])] = 0
      B[numpy.random.randint(B.shape[0]),numpy.random.randint(B.shape[1])] = 0

    A = ca.sparsify(A)
    B = ca.sparsify(B)

    self.checkarray(ca.norm_inf_mul(A,B),ca.norm_inf(A @ B))
    self.checkarray(ca.DM(ca.norm_0_mul(A,B)),mtimes(A,B).nnz())


    A = numpy.random.random((8,2))
    B = numpy.random.random((2,10))

    self.checkarray(ca.norm_inf_mul(A,B),ca.norm_inf(A @ B))
    self.checkarray(ca.DM(ca.norm_0_mul(A,B)),mtimes(A,B).nnz())

    # Sparse
    for i in range(5):
      A[numpy.random.randint(A.shape[0]),numpy.random.randint(A.shape[1])] = 0
      B[numpy.random.randint(B.shape[0]),numpy.random.randint(B.shape[1])] = 0

    A = ca.sparsify(A)
    B = ca.sparsify(B)

    self.checkarray(ca.norm_inf_mul(A,B),ca.norm_inf(A @ B))
    self.checkarray(ca.DM(ca.norm_0_mul(A,B)),mtimes(A,B).nnz())

  def  test_mul3_issue_1465(self):

    w2 = ca.SX.sym("w",2,1)
    w3 = ca.SX.sym("w",3,1)
    Q = np.eye(2)

    mtimes(w2.T,Q)
    w2.T @ Q
    with self.assertInException("Matrix product with incompatible dimensions. Lhs is 1x3 and rhs is 2x2."):
      mtimes(w3.T,Q)
    with self.assertInException("Matrix product with incompatible dimensions. Lhs is 1x3 and rhs is 2x2."):
      w3.T @ Q

    mtimes([w2.T,Q,w2])
    with self.assertInException("Matrix product with incompatible dimensions. Lhs is 1x3 and rhs is 2x2."):
      mtimes([w3.T,Q,w2])
    with self.assertInException("Matrix product with incompatible dimensions. Lhs is 1x2 and rhs is 3x1."): 
      mtimes([w2.T,Q,w3])
    with self.assertInException("Matrix product with incompatible dimensions. Lhs is 1x3 and rhs is 2x2."):
      mtimes([w3.T,Q,w3])

    with self.assertInException("unsupported operand type(s) for @: 'SX' and 'NoneType'"):
      w2.T @ None # pyright: ignore[reportCallIssue,reportArgumentType, reportOperatorIssue]
    
    with self.assertInException("mtimes(Sparsity,Sparsity,str)"):
      mtimes(w2,None) # pyright: ignore[reportCallIssue,reportArgumentType]
  
    with self.assertInException("mtimes(Sparsity,Sparsity,str)"):
      mtimes(w2.T,Q,None) # pyright: ignore[reportCallIssue,reportArgumentType]

    with self.assertInException("mtimes(Sparsity,Sparsity,str)"):
      mtimes(w2.T,Q,w2) # pyright: ignore[reportCallIssue,reportArgumentType]

  def test_chol(self):
    numpy.random.seed(0)

    for i in range(4):
      A = numpy.random.random((3,3))
      H = A @ A.T

      R = ca.chol(H)

      assert R.is_triu()
      self.checkarray(R.T @ R,H)
  def test_skew(self):
    x = ca.DM([1,7,13])
    self.checkarray(ca.inv_skew(ca.skew(x)),x)
    y = ca.DM([0.2,0.9,0.4])
    self.checkarray(ca.skew(x) @ y,ca.cross(x,y))

  def test_nz_overflow(self):
    d = ca.DM([2,3])
    r = d.nz[:]
    self.checkarray(r,d)

  def test_DMcrash(self):
    with self.assertRaises(Exception):
      ca.DM([ca.DM([1,2]),ca.DM([1,2])])  # pyright: ignore[reportCallIssue,reportArgumentType]
    a = ca.DM([ca.DM([1]),ca.DM([2])])  # pyright: ignore[reportCallIssue,reportArgumentType]
    self.checkarray(a,ca.DM([1,2]))

  def test_sparsity_operation(self):
    L = [ca.DM(1), ca.DM(ca.Sparsity(1,1),1), ca.DM(ca.Sparsity(2,1),1), ca.DM(ca.Sparsity.dense(2,1),1)]

    for a in L:
      for b in L:
        c = a*b

        if a.nnz()==0 or b.nnz()==0:
          self.assertTrue(c.nnz()==0)
        else:
          self.assertTrue(c.nnz()>0)

    self.assertTrue(ca.sum2(ca.DM(ca.Sparsity(1,1),1)).nnz()==0)

  def test_matlab_operations(self):

    data = [ np.array([[1,3],[11,17]]) , np.array([[1,3]]) ,np.array([[1],[3]]), np.array([[3]])]

    for A in data:
      B = ca.reshape(ca.DM(A),A.shape)  # pyright: ignore[reportCallIssue,reportArgumentType]
      #self.checkarray(np.cumsum(A),cumsum(B))
      self.checkarray(np.cumsum(A,0),ca.cumsum(B,0))
      self.checkarray(np.cumsum(A,1),ca.cumsum(B,1))

      Bs = ca.MX.sym("B",B.shape)
      A0 = ca.cumsum(B,0)
      A1 = ca.cumsum(B,1)

      self.checkarray(np.cumsum(A,0),ca.evalf(ca.cumsum(ca.MX(B),0)))
      self.checkarray(np.cumsum(A,1),ca.evalf(ca.cumsum(B,1)))

      f = ca.Function("f",[Bs],[A0,A1])
      self.checkfunction(f,f.expand(),inputs = [B])
      self.check_codegen(f,inputs=[B],check_serialize=True)

      #self.checkarray(np.diff(A),diff(B))

      #self.checkarray(np.diff(A,1),diff(B,1))
      #self.checkarray(np.diff(A,2),diff(B,2))

      self.checkarray(np.diff(A,1,0),ca.diff(B,1,0))
      #self.checkarray(np.diff(A,1,1),diff(B,1,1))


  def test_singular_repmat(self):
    for X in [ca.DM, ca.SX, ca.MX, ca.Sparsity]:
      for n_b in [0,2]:
        for m_b in [0,2]:
          b = X(n_b, m_b)

          for n in [0,3]:
            for m in [0,3]:
              self.assertEqual(ca.repmat(b, n, m).shape,(n_b * n, m_b * m))


  def test_serialize(self):
    for a in [ca.DM(), ca.DM(2), ca.DM([1,2]),ca.DM([[1,2],[3,4],[5,6]])]:
      b = ca.DM.deserialize(a.serialize())
      self.checkarray(a,b)

  def test_iterable(self):
    a = ca.DM([1,2,3])
    b = list(iter(a.nz))
    self.checkarray(a,ca.DM(b))  # pyright: ignore[reportCallIssue,reportArgumentType]

    with self.assertInException("CasADi matrices are not iterable"):
      iter(a)

  def test_from_file(self):
    a = ca.sparsify(ca.DM([[1,0,-6,5,0],[4,0,-4e-301,9.3e-18,0]]))

    a.print_dense()

    a[0,4] = 0
    a[1,1] = 0
    a.to_file("test.mtx")

    with self.assertInException("foobar.mtx"):
      ca.DM.from_file("foobar.mtx")
    b = ca.DM.from_file("test.mtx")

    self.assertTrue(a.sparsity()==b.sparsity())
    self.checkarray(a,b)

    a.to_file("test.txt")
    b = ca.DM.from_file("test.txt")
    self.checkarray(a,b)

    for s in ["1 00 -6 5 0\n4 0 -4e-301 9.3e-18 00\n",
              "1 00 -6 5 0\n4 0 -4e-301 9.3e-18 00",
              "1 00 -6 5 0 \n4 0 -4e-301 9.3e-18 00   ",
              "1\t00\t-6\t5\t0\t\n4 0 -4e-301 9.3e-18 00\n",
              "1 00 -6 5 0\n%foobar\n4 0 -4e-301 9.3e-18 00\n"]:

      with open("test2.txt","w") as f:
        f.write(s)
      b = ca.DM.from_file("test2.txt")
      self.assertTrue(a.sparsity()==b.sparsity())
      self.assertTrue(float(ca.norm_1(a-b))==0)

    with open("test.txt","w") as f:
      f.write("1 aa 6\n4 8 9e-18\n")


    with self.assertInException("Parsing"):
      ca.DM.from_file("test.txt")

    a = ca.DM([3, inf, -inf, np.nan])
    a.to_file("test.txt")
    b = ca.DM.from_file("test.txt")
    self.checkarray(a,b)
  def test_norm_2(self):
    a = np.array([1,2,3])
    r = np.linalg.norm(a)
    for M in [ca.DM,ca.MX,ca.SX]:
      self.checkarray(ca.evalf(ca.norm_2(M(a))),r)
      self.checkarray(ca.evalf(ca.norm_2(M(a).T)),r)

  def test_ldl(self):
    numpy.random.seed(1)
    ca.DM.rng(1)
    H = ca.diagcat(ca.DM.rand(5,5),ca.DM.rand(5,5))
    H = H+H.T+2*ca.DM.eye(10)

    p = np.random.permutation(10)
    H = H[p,p]


    [D,Lt,p] = ca.ldl(H)
    P = ca.DM.eye(10)[:,p]

    F = ca.sqrt(ca.diag(D)) @ (ca.DM.eye(10)+Lt) @ P.T
    print(H-F.T @ F)
    self.assertTrue(ca.norm_fro(H-F.T @ F)<=1e-14)


  def test_im_bugs(self):
    a = ca.vertcat(1,2)
    self.assertTrue(isinstance(a,ca.DM))
    self.checkarray(c.linspace(1,3,10),c.linspace(1.0,3.0,10))

  def test_linspace(self):
    # CasADi's linspace must agree bit-for-bit with numpy.linspace, for
    # DM (numeric) and for MX/SX after Function evaluation. This pins down
    # the FP recipe so future refactors (incl. MX node-count optimizations)
    # cannot silently change the numerics.
    cases = [(1.0, 3.0), (-2.5, 7.25), (0.0, 1.0), (1e10, 1e-10), (-1.0, 1.0)]
    sizes = [2, 3, 5, 10, 100, 1001]
    for a, b in cases:
      for n in sizes:
        ref = numpy.linspace(a, b, n)

        dm = numpy.array(c.linspace(ca.DM(a), ca.DM(b), n)).ravel()
        self.assertTrue(numpy.array_equal(dm, ref),
                        "DM linspace mismatch a=%g b=%g n=%d" % (a, b, n))

        x = ca.MX.sym('x'); y = ca.MX.sym('y')
        f = ca.Function('f', [x, y], [c.linspace(x, y, n)])
        mx = numpy.array(f(a, b)).ravel()
        self.assertTrue(numpy.array_equal(mx, ref),
                        "MX linspace mismatch a=%g b=%g n=%d" % (a, b, n))

        x = ca.SX.sym('x'); y = ca.SX.sym('y')
        f = ca.Function('f', [x, y], [c.linspace(x, y, n)])
        sx = numpy.array(f(a, b)).ravel()
        self.assertTrue(numpy.array_equal(sx, ref),
                        "SX linspace mismatch a=%g b=%g n=%d" % (a, b, n))

  def test_linspace_mx_node_count(self):
    # MX linspace should yield O(1) graph nodes regardless of nsteps.
    x = ca.MX.sym('x'); y = ca.MX.sym('y')
    n_small = ca.n_nodes(c.linspace(x, y, 5))
    n_large = ca.n_nodes(c.linspace(x, y, 1000))
    # Allow some slack but reject O(n) growth.
    self.assertTrue(n_large < n_small + 5,
                    "MX linspace node count grew with nsteps: %d -> %d"
                    % (n_small, n_large))

  def test_permutation(self):
    n = 10
    numpy.random.seed(1)
    p = np.random.permutation(n)
    S = ca.Sparsity.permutation(p)
    self.checkarray(ca.DM(p).T,S.permutation_vector())
    self.assertTrue(S.is_permutation())
    v = ca.DM.rand(n)
    self.checkarray(S @ v, v[p])
    S = ca.Sparsity.permutation(p, True)
    self.checkarray(((S @ v))[p], v)

  def test_sparsity_orthonormality(self):
    alltests = [('is_orthonormal_rows',True),('is_orthonormal_columns',True),('is_selection',True),('is_orthonormal_rows',False),('is_orthonormal_columns',False),('is_selection',False),('is_permutation',),('is_orthonormal',False),('is_orthonormal',True)]
    
    def check_properties(S,prop):
      for p in prop:
        assert p in alltests
      for t in alltests:
        quality = t in prop
        method = getattr(S,t[0])
        result = method(*t[1:])
        self.assertEqual(quality,result,"Call %s expected to return %s" % (str(t),str(quality)))
        
    for S in [ca.sparsify(ca.DM([[0,0,0,0,0],
                           [0,0,0,0,0],
                           [0,0,0,0,0],
                           [0,0,0,0,0]])),
              ca.sparsify(ca.DM([[0,0,0,0,0],
                           [0,1,0,0,0],
                           [0,0,0,0,0],
                           [0,0,0,0,0]])),
              ca.sparsify(ca.DM([[0,0,0,0,0],
                           [0,1,0,0,0],
                           [0,0,0,1,0],
                           [0,0,0,0,0]]))]:
      check_properties(S.sparsity(), [('is_orthonormal_rows',True),('is_orthonormal_columns',True),('is_selection',True),('is_orthonormal',True)]) 

    for S in [ca.sparsify(ca.DM([[0,0,0,0,0],
                           [0,1,1,0,0],
                           [0,0,0,0,0],
                           [0,0,0,0,0]])),
              ca.sparsify(ca.DM([[0,1,0,0],
                           [1,0,0,0],
                           [0,1,1,0],
                           [0,0,0,0],
                           [0,0,0,1]])),
              ca.sparsify(ca.DM([[0,0,0,0,0],
                           [0,1,0,0,0],
                           [0,1,0,0,0],
                           [0,0,0,0,0]]))]:
      check_properties(S.sparsity(), []) 
    
    S = ca.sparsify(ca.DM([[0,1,0,0,0],
                     [1,0,0,0,0],
                     [0,0,1,0,0],
                     [0,0,0,1,0]]))
    check_properties(S.sparsity(), [('is_orthonormal_rows',True),('is_orthonormal_columns',True),('is_selection',True),('is_orthonormal_rows',False),('is_selection',False),('is_orthonormal',True)])
    
    S = ca.sparsify(ca.DM([[0,1,0,0],
                     [1,0,0,0],
                     [0,0,1,0],
                     [0,0,0,1]]))
    check_properties(S.sparsity(), [('is_orthonormal_rows',True),('is_orthonormal_columns',True),('is_selection',True),('is_orthonormal_rows',False),('is_orthonormal_columns',False),('is_selection',False),('is_permutation',),('is_orthonormal',True),('is_orthonormal',False)])

    S = ca.sparsify(ca.DM([[0,1,0,0],
                     [1,0,0,0],
                     [0,0,1,0],
                     [0,0,0,0],
                     [0,0,0,1]]))
    check_properties(S.sparsity(), [('is_orthonormal_rows',True),('is_orthonormal_columns',True),('is_selection',True),('is_orthonormal_columns',False),('is_orthonormal',True)])
    
  def test_permutation_solve(self):
    n = 10
    numpy.random.seed(1)
    ca.DM.rng(1)
    p = np.random.permutation(n)
    P = ca.Sparsity.permutation(p)
    fs = []
    for X in [ca.SX,ca.MX]:
      x = X.sym("x",P)
      b = X.sym("b",n)
      assert x.sparsity().is_orthonormal()
      nz = ca.sparsity_cast(x, ca.Sparsity.dense(n))
      f = ca.Function("f",[x,b],[ca.solve(x,b)])
      f.generate('f.c')
      fs.append(f)    
    x0 = ca.DM.rand(x.sparsity())
    b0 = ca.DM.rand(b.sparsity())
    res0 = fs[0](x0,b0)
    res1 = fs[1](x0,b0)
    self.checkarray(res0,res1)

  def test_permutation_inv(self):
    n = 10
    numpy.random.seed(1)
    ca.DM.rng(1)
    p = np.random.permutation(n)
    P = ca.Sparsity.permutation(p)
    fs = []
    for X in [ca.SX,ca.MX]:
      x = X.sym("x",P)
      f = ca.Function("f",[x],[ca.solve(x,X.eye(n))])
      f.generate('f.c')
      fs.append(f)    
    x0 = ca.DM(P,ca.DM.rand(n))
    res0 = fs[0](x0)
    res1 = fs[1](x0)
    self.checkarray(res0,res1)

    fs = []
    for X in [ca.SX,ca.MX]:
      x = X.sym("x",P)
      f = ca.Function("f",[x],[ca.inv(x)])
      f.disp(True)
      fs.append(f)    
    x0 = ca.DM(P,ca.DM.rand(n))
    res0 = fs[0](x0)
    res1 = fs[1](x0)
    self.checkarray(res0,res1)

    
  def test_mac_simplify(self):
    n = 10
    numpy.random.seed(1)
    ca.DM.rng(1)
    p = np.random.permutation(n)
    P = ca.Sparsity.permutation(p)
    fs = []
    for X in [ca.SX,ca.MX]:
      b = X.sym("b",n)
      x = X.sym("x",P)
      f = ca.Function("f",[x,b],[x @ b])
      f.generate('f.c')
      fs.append(f)  
    x0 = ca.DM(P,ca.DM.rand(n))
    b0 = ca.DM.rand(b.sparsity())
    res0 = fs[0](x0,b0)
    res1 = fs[1](x0,b0)
    self.checkarray(res0,res1)
    

  def test_inplace(self):
    n = 10
    x = ca.MX.sym("x",n)
    y = ca.MX(ca.Sparsity.diag(n),x)
    f = ca.Function('f',[x],[y])
    self.assertTrue(f.sz_w()==n)


  def test_horzsplit_n(self):
    self.assertTrue(len(ca.horzsplit_n(ca.DM.rand(1,6),2))==2)
    self.assertTrue(len(ca.vertsplit_n(ca.DM.rand(6),2))==2)
    self.assertTrue(len(ca.horzsplit_n(ca.DM(0,0),2))==2)
    self.assertTrue(len(ca.vertsplit_n(ca.DM(0,0),2))==2)
    with self.assertRaises(Exception):
       ca.horzsplit_n(ca.DM.rand(1,5),2)
    with self.assertRaises(Exception):
       ca.vertsplit_n(ca.DM.rand(5),2)
       
if __name__ == '__main__':
    unittest.main()
