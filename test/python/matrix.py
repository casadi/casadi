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
    a=DM(array([[1,2,3],[4,5,6],[7,8,9]]))
    b=DM([[1,2,3],[4,5,6],[7,8,9]])
    self.checkarray(a,b,"List of list constructor")

  def test_sum(self):
    self.message("sum")
    D=DM([[1,2,3],[4,5,6],[7,8,9]])
    self.checkarray(c.sum1(D),array([[12,15,18]]),'sum()')
    self.checkarray(c.sum2(D),array([[6,15,24]]).T,'sum()')


  def test_inv(self):
    self.message("Matrix inverse")
    a = DM([[1,2],[1,3]])
    self.checkarray(mtimes(c.inv(a),a),eye(2),"DM inverse")

  def test_trans(self):
    self.message("trans")
    a = DM(0,1)
    b = a.T
    self.assertEqual(b.size1(),1)
    self.assertEqual(b.size2(),0)

  def test_vertcat(self):
    self.message("vertcat")
    A = DM.ones(2,3)
    B = DM(4,3)
    C = vertcat(A,B)

    self.checkarray(C.shape,(6,3),"vertcat shape")
    self.assertEqual(C.nnz(),A.nnz(),"vertcat size")

    self.assertRaises(RuntimeError,lambda : horzcat(A,B))

  def test_horzcat(self):
    self.message("horcat")
    A = DM.ones(3,2)
    B = DM(3,4)
    C = horzcat(A,B)

    self.checkarray(C.shape,(3,6),"horzcat shape")
    self.assertEqual(C.nnz(),A.nnz(),"vertcat size")

    self.assertRaises(RuntimeError,lambda : vertcat(A,B))

  def test_diagcat(self):

    x = MX.sym("x",2,2)
    y = MX.sym("y",Sparsity.lower(3))
    z = MX.sym("z",4,2)

    L = [x,y,z]

    fMX = Function("fMX", L,[diagcat(*L)])

    LSX = [ SX.sym("",i.sparsity()) for i in L ]
    fSX = Function("fSX", LSX,[diagcat(*LSX)])

    for f in [fMX,fSX]:
      for i in range(3):
        f_in[i]=list(range(f.nnz_in(i)))

    self.checkfunction(fMX,fSX)

  def test_veccat(self):
    self.message("vecccat")
    A = DM(2,3)
    A[0,1] = 2
    A[1,0] = 1
    A[1,2] = 3
    B = DM(3,1)
    B[0,0] = 4
    B[1,0] = 5
    B[2,0] = 6
    C = veccat(*[A,B])

    self.checkarray(C.shape,(9,1),"veccat shape")
    self.assertEqual(C.nnz(),A.nnz()+B.nnz(),"veccat size")

    self.checkarray(tuple(C.nonzeros()),tuple(arange(1,7)),"numbers shape")

  def test_slicestepnegative(self):
    self.message("Slice step negative")
    a1 = [1,2,3,4,5]
    a2 = DM(a1)

    self.checkarray(a2[0:4:-1,0],DM(a1[0:4:-1])) # gives empty set
    self.checkarray(a2[4:0:-1,0],DM(a1[4:0:-1])) # gives [5, 4, 3, 2]
    self.checkarray(a2[0:4:-2,0],DM(a1[0:4:-2])) # gives empty set
    self.checkarray(a2[4:0:-2,0],DM(a1[4:0:-2])) # gives [5, 4, 3, 2]
    self.checkarray(a2[1:4:-2,0],DM(a1[1:4:-2])) # gives empty set
    self.checkarray(a2[4:1:-2,0],DM(a1[4:1:-2])) # gives [5, 4, 3, 2]
    self.checkarray(a2[0:3:-2,0],DM(a1[0:3:-2])) # gives empty set
    self.checkarray(a2[3:0:-2,0],DM(a1[3:0:-2])) # gives [5, 4, 3, 2]
    self.checkarray(a2[::-1,0],DM(a1[::-1])) # gives [5, 4, 3, 2, 1]
    self.checkarray(a2[::1,0],DM(a1[::1])) # gives [1,2,3,4,5]
    self.checkarray(a2[2::-1,0],DM(a1[2::-1])) # gives [3,2,1]
    self.checkarray(a2[:2:-1,0],DM(a1[:2:-1])) # gives [5,4]
    self.checkarray(a2[-1::-1,0],DM(a1[-1::-1])) # gives [5,4,3,2,1]
    self.checkarray(a2[:-1:-1,0],DM(a1[:-1:-1])) # gives []

  def test_indexingOutOfBounds(self):
    self.message("Indexing out of bounds")
    y = DM.zeros(4, 5)
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
    a = SX.sym("a",Sparsity.diag(50000))

    a[:,:]

  def test_nonmonotonous_indexing(self):
    self.message("non-monotonous indexing")
    # Regression test for #354
    A = DM([[1,2,3],[4,5,6],[7,8,9]])
    B = A[[0,2,1],0]
    self.checkarray(DM([1,7,4]),B,"non-monotonous")

    B = A[0,[0,2,1]]
    self.checkarray(DM([1,3,2]).T,B,"non-monotonous")

  def test_IM_indexing(self):
    self.message("IM")
    A = IM(2,2)
    A[0,0] = 1
    A[1,1] = 3
    A[0,1] = 2
    A[1,0] = 4


    B = DM([1,2,3,4,5])

    B_ = B[A]

    self.checkarray(B_,DM([[2,3],[5,4]]),"Imatrix indexing")

    B[A] = DM([[1,2],[3,4]])

    self.checkarray(B,DM([1,1,2,4,3]),"Imatrix indexing assignement")

    #B[A].set(DM([[10,20],[30,40]]))

    #self.checkarray(B,DM([1,10,20,40,30]),"Imatrix indexing setting")

    B = IM([1,2,3,4,5])

    B_ = B[A]

    self.checkarray(array(B_),DM([[2,3],[5,4]]),"Imatrix indexing")

    B[A] = IM([[1,2],[3,4]])

    self.checkarray(array(B),DM([1,1,2,4,3]),"Imatrix indexing assignement")

    B = DM(5,1)

    self.assertRaises(Exception, lambda : B.nz[A])

  def test_sparsity_indexing(self):
    self.message("sparsity")

    B = DM([[1,2,3,4,5],[6,7,8,9,10]])

    A = IM([[1,1,0,0,0],[0,0,1,0,0]])
    A = sparsify(A)
    sp = A.sparsity()


    self.checkarray(B[sp],DM([[1,2,0,0,0],[0,0,8,0,0]]),"sparsity indexing")

    B[sp] = -4

    self.checkarray(B,DM([[-4,-4,3,4,5],[6,7,-4,9,10]]),"sparsity indexing assignement")

    B = DM([[1,2,3,4,5],[6,7,8,9,10]])

    B[sp] = 2*B

    self.checkarray(B,DM([[2,4,3,4,5],[6,7,16,9,10]]),"Imatrix indexing assignement")

    self.assertRaises(Exception, lambda : B[Sparsity.dense(4,4)])

  def test_index_setting(self):
    self.message("index setting")
    B = DM([1,2,3,4,5])

    B[0] = 8
    self.checkarray(B,DM([8,2,3,4,5]),"index setting")
    B[1,0] = 4
    self.checkarray(B,DM([8,4,3,4,5]),"index setting")
    B[:,0] = 7
    self.checkarray(B,DM([7,7,7,7,7]),"index setting")
    #B[0].set(3)
    #self.checkarray(B,DM([3,7,7,7,7]),"index setting")
    #B[0].setAll(4)
    #self.checkarray(B,DM([4,4,4,4,4]),"index setting")


  def test_issue298(self):
    self.message("Issue #298")
    a = DM(4,1)
    b = c.reshape(a,2,2)
    self.assertEqual(type(a),type(b))

    a = IM(4,1)
    b = DM(a)
    self.assertTrue(isinstance(b,DM))

    a = DM(4,1)
    b = IM(a)
    self.assertTrue(isinstance(b,IM))

  def test_det(self):
    self.message("Determinant")
    npy_det = numpy.linalg.det

    a = DM(1,1)
    a[0,0] = 5
    self.checkarray(det(a),npy_det(a),"det()")

    a = DM(5,5)
    for i in range(5):
      a[i,i] = i+1

    self.checkarray(det(a),npy_det(a),"det()")

    a = DM(5,5)
    for i in range(4):
      a[i,i] = i+1
    a[0,4] = 3
    a[4,0] = 7

    self.checkarray(det(a),npy_det(a),"det()")

    a = DM(5,5)
    for i in range(5):
      for j in range(5):
        a[i,j] = i+j

    self.checkarray(det(a),npy_det(a),"det()")

    a = DM(5,5)
    for i in range(4):
      for j in range(5):
        a[i,j] = i+j

    self.checkarray(det(a),npy_det(a),"det()")

    a = DM(5,5)
    for i in range(5):
      for j in range(4):
        a[i,j] = i+j

    self.checkarray(det(a),npy_det(a),"det()")

    a = DM(5,5)
    for i in range(4):
      for j in range(5):
        a[i,j] = i+j
    a[4,1] = 12

    self.checkarray(det(a),npy_det(a),"det()")

    a = DM(5,5)
    for i in range(5):
      for j in range(4):
        a[i,j] = i+j
    a[1,4] = 12

    self.checkarray(det(a),npy_det(a),"det()")

    a = DM(5,5)
    for i in range(4):
      for j in range(5):
        a[i,j] = i+j
    a[4,2] = 12

    self.checkarray(det(a),npy_det(a),"det()")

    a = DM(5,5)
    for i in range(5):
      for j in range(4):
        a[i,j] = i+j
    a[2,4] = 12

    self.checkarray(det(a),npy_det(a),"det()")

    a = DM(50,50)
    for i in range(50):
      a[i,i] = i+1

    self.checkarray(det(a)/npy_det(a),1,"det()")

  @skip(not GlobalOptions.getSimplificationOnTheFly())
  def test_inv_sparsity(self):
    self.message("sparsity pattern of inverse")

    n = 8

    sp = Sparsity.lower(n)

    x  = SX(sp,vertcat(*[SX.sym("a%d" % i) for i in range(sp.nnz())]))


    x_ = DM.ones(x.sparsity())

    I_ = DM.ones(inv(x).sparsity())

    # For a reducible matrix, struct(A^(-1)) = struct(A)
    self.checkarray(x_,I_,"inv")

    sp = Sparsity.lower(n)

    x = SX.sym("a", sp)
    x[0,n-1] = 1


    I_ = DM.ones(inv(x).sparsity())

    # An irreducible matrix has a dense inverse in general
    self.checkarray(DM.ones(n,n),I_,"inv")

    x = SX.sym("a", sp)
    x[0,int(n/2)] = 1

    s_ = DM.ones(sp)
    s_[:,:int(n/2)+1] = 1

    I_ = DM.ones(inv_minor(x).sparsity())

    s_ = densify(s_)
    T_ = densify(I_)
    # An irreducible matrix does not have to be dense per se
    self.checkarray(s_,I_,"inv")

  def test_Imatrix_operations(self):
    self.message("IM operations")
    a = IM.ones(2,2)
    b = horzcat(a,a)
    self.assertTrue(isinstance(b,IM))

  def test_mtimes(self):
    A = DM.ones((4,3))
    B = DM.ones((3,8))
    C = DM.ones((8,7))

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
    B = DM([[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16],[17,18,19,20]])

    A = DM(B)
    A.remove([],[])
    self.checkarray(A, B,"remove nothing")

    A = DM(B)
    A.remove([],[1])
    self.checkarray(A, DM([[1,3,4],[5,7,8],[9,11,12],[13,15,16],[17,19,20]]),"remove a column")

    A = DM(B)
    A.remove([0,3],[1])
    self.checkarray(A, DM([[5,7,8],[9,11,12],[17,19,20]]),"remove a column and two rows ")

  def test_comparisons(self):
    for m in [DM,IM]:
      A = m([[5,4],[2,1]])

      for c in [6,6.0,DM([6]),IM([6]),matrix(6)]:
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
        if args.known_bugs or not isinstance(c,matrix):
          self.checkarray(c==A,m([[0,0],[0,0]]),"==")
          self.checkarray(c!=A,m([[1,1],[1,1]]),"!=")

      for c in [5,5.0,DM([5]),IM([5]),matrix(5)]:
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
        if args.known_bugs or not isinstance(c,matrix):
          self.checkarray(c==A,m([[1,0],[0,0]]),"==")
          self.checkarray(c!=A,m([[0,1],[1,1]]),"!=")

      for c in [4,4.0,DM([4]),IM([4]),matrix(4)]:
        self.checkarray(A<=c,m([[0,1],[1,1]]),"<=")
        self.checkarray(A<c,m([[0,0],[1,1]]),"<")
        self.checkarray(A>c,m([[1,0],[0,0]]),">")
        self.checkarray(A>=c,m([[1,1],[0,0]]),">=")
        if args.known_bugs or not isinstance(c,matrix):
          self.checkarray(A==c,m([[0,1],[0,0]]),"==")
          self.checkarray(A!=c,m([[1,0],[1,1]]),"!=")

        self.checkarray(c>=A,m([[0,1],[1,1]]),"<=")
        self.checkarray(c>A,m([[0,0],[1,1]]),"<")
        self.checkarray(c<A,m([[1,0],[0,0]]),">")
        self.checkarray(c<=A,m([[1,1],[0,0]]),">=")
        if args.known_bugs or not isinstance(c,matrix):
          self.checkarray(c==A,m([[0,1],[0,0]]),"==")
          self.checkarray(c!=A,m([[1,0],[1,1]]),"!=")

      for c in [1,1.0,DM([1]),IM([1]),matrix(1)]:
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
        if args.known_bugs or not isinstance(c,matrix):
          self.checkarray(c==A,m([[0,0],[0,1]]),"==")
          self.checkarray(c!=A,m([[1,1],[1,0]]),"!=")

      for c in [0,DM([0]),IM([0]),matrix(0)]:
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
        if args.known_bugs or not isinstance(c,matrix):
          self.checkarray(c==A,m([[0,0],[0,0]]),"==")
          self.checkarray(c!=A,m([[1,1],[1,1]]),"!=")

  def test_truth(self):
    self.assertTrue(bool(DM([1])))
    self.assertFalse(bool(DM([0])))
    self.assertTrue(bool(DM([0.2])))
    self.assertTrue(bool(DM([-0.2])))
    self.assertRaises(Exception, lambda : bool(DM([2.0,3])))
    self.assertRaises(Exception, lambda : bool(DM()))

  def test_listslice(self):
    def check(d,rowbase,colbase):
      for col in permutations(colbase):
        for row in permutations(rowbase):
          r = IM.zeros(len(row),len(col))
          for i,ii in enumerate(row):
            for j,jj in enumerate(col):
              r[i,j] = d[ii,jj]
          self.checkarray(d[row,col],r,"%s[%s,%s]" % (repr(d),str(row),str(col)))


    # get1
    check(IM(Sparsity.dense(3,3),list(range(3*3))),[0,1,2],[0,1,2])
    check(IM(Sparsity.dense(4,4),list(range(4*4))),[0,1,3],[0,2,3])
    check(IM(Sparsity.dense(3,3),list(range(3*3))),[0,0,1],[0,0,1])
    check(IM(Sparsity.dense(3,3),list(range(3*3))),[0,0,2],[0,0,2])
    check(IM(Sparsity.dense(3,3),list(range(3*3))),[1,1,2],[1,1,2])

    sp = Sparsity.lower(4)
    d = IM(sp,list(range(sp.nnz())))
    check(d,[0,1,3],[0,2,3])
    check(d.T,[0,1,3],[0,2,3])

    sp = Sparsity.rowcol([0,1,2],[0,1],4,4)
    d = IM(sp,list(range(sp.nnz())))
    check(d,[0,3],[0,2])

    # get2
    check(IM(Sparsity.dense(2,2),list(range(2*2))),[0,0,0],[0,0,0])
    check(IM(Sparsity.dense(2,2),list(range(2*2))),[0,0,1],[0,0,1])
    check(IM(Sparsity.dense(2,2),list(range(2*2))),[1,1,0],[1,1,0])
    check(IM(Sparsity.dense(2,2),list(range(2*2))),[1,1,1],[1,1,1])

    sp = Sparsity.lower(3)
    d = IM(sp,list(range(sp.nnz())))
    check(d,[0,1,2],[0,1,2])
    check(d.T,[0,1,2],[0,1,2])

    sp = Sparsity.rowcol([0,2],[0,1],4,4)
    d = IM(sp,list(range(sp.nnz())))
    check(d,[0,1,3],[0,2,3])

  def test_sparsesym(self):
    # feature removed in 73f407e
    return
    self.message("sparsesym")
    D = DM([[1,2,-3],[2,-1,0],[-3,0,5]])
    D = sparsify(D)

    i = D.getSym()
    #self.checkarray(list(i),[1,2,-1,-3,5])
    A = 2*D
    A.setSym(i)
    self.checkarray(A,D)

  def test_diagcat(self):
    self.message("diagcat")
    C = diagcat(*[DM([[-1.4,-3.2],[-3.2,-28]]),DM([[15,-12,2.1],[-12,16,-3.8],[2.1,-3.8,15]]),1.8,-4.0])
    r = DM([[-1.4,-3.2,0,0,0,0,0],[-3.2,-28,0,0,0,0,0],[0,0,15,-12,2.1,0,0],[0,0,-12,16,-3.8,0,0],[0,0,2.1,-3.8,15,0,0],[0,0,0,0,0,1.8,0],[0,0,0,0,0,0,-4]])
    r = sparsify(r)
    self.checkarray(C,r)

  def test_diag_sparse(self):
    self.message("diag sparse")

    for n in [[0,1,0,0,2,3,4,5,6,0],[1,2,3,0],[0,1,2,3]]:
      d = DM(n)
      D = DM(n)
      d = sparsify(d)
      m = c.diag(d)
      M = sparsify(c.diag(D))

      self.checkarray(m.sparsity().colind(),M.sparsity().colind())
      self.checkarray(m.sparsity().row(),M.sparsity().row())

  def test_sprank(self):
    self.message("sprank")

    a = DM([[1,0,0],[0,1,0],[0,0,1]])
    a = sparsify(a)
    self.assertEqual(sprank(a),3)

    a = DM([[1,0,0],[0,0,0],[0,0,1]])
    a = sparsify(a)
    self.assertEqual(sprank(a),2)

    a = DM([[0,0,0],[0,0,0],[0,0,1]])
    a = sparsify(a)
    self.assertEqual(sprank(a),1)

    a = DM([[0,0,0],[0,0,0],[0,0,0]])
    a = sparsify(a)
    self.assertEqual(sprank(a),0)

    self.assertEqual(sprank(DM.ones(1,3)),1)
    self.assertEqual(sprank(DM.ones(3,1)),1)
    self.assertEqual(sprank(DM.ones(2,3)),2)
    self.assertEqual(sprank(DM.ones(3,2)),2)
    self.assertEqual(sprank(DM.ones(3,3)),3)
    self.assertEqual(sprank(DM.ones(3,3)),3)

    A = DM(6,4)
    A[0,0] = 1
    A[1,2] = 1
    A[2,2] = 1
    A[5,3] = 1

    self.assertEqual(sprank(A),3)

  def test_cross(self):
    self.message("cross products")

    crossc = c.cross

    self.checkarray(crossc(DM([1,0,0]),DM([0,1,0])),DM([0,0,1]))

    self.checkarray(crossc(DM([1.1,1.3,1.7]),DM([2,3,13])),DM([11.8,-10.9,0.7]))
    self.checkarray(crossc(DM([1.1,1.3,1.7]).T,DM([2,3,13]).T),DM([11.8,-10.9,0.7]).T)

    self.checkarray(crossc(DM([[1.1,1.3,1.7],[1,0,0],[0,0,1],[4,5,6]]),DM([[2,3,13],[0,1,0],[0,0,1],[1,0,1]])),DM([[11.8,-10.9,0.7],[0,0,1],[0,0,0],[5,2,-5]]))
    self.checkarray(crossc(DM([[1.1,1.3,1.7],[1,0,0],[0,0,1],[4,5,6]]).T,DM([[2,3,13],[0,1,0],[0,0,1],[1,0,1]]).T),DM([[11.8,-10.9,0.7],[0,0,1],[0,0,0],[5,2,-5]]).T)

    self.checkarray(crossc(DM([[1.1,1.3,1.7],[1,0,0],[0,0,1],[4,5,6]]),DM([[2,3,13],[0,1,0],[0,0,1],[1,0,1]]),2),DM([[11.8,-10.9,0.7],[0,0,1],[0,0,0],[5,2,-5]]))

    self.checkarray(crossc(DM([[1.1,1.3,1.7],[1,0,0],[0,0,1],[4,5,6]]).T,DM([[2,3,13],[0,1,0],[0,0,1],[1,0,1]]).T,1),DM([[11.8,-10.9,0.7],[0,0,1],[0,0,0],[5,2,-5]]).T)

  def test_is_regular(self):
    self.assertTrue(DM([1,2]).is_regular())
    self.assertFalse(DM([1,inf]).is_regular())
    self.assertFalse(DM.nan(2).is_regular())

  def test_sizes(self):
    self.assertEqual(Sparsity.diag(10).nnz_diag(),10)
    self.assertEqual(Sparsity.diag(10).nnz_upper(),10)
    self.assertEqual(Sparsity.diag(10).nnz_lower(),10)
    self.assertEqual(Sparsity.dense(10,10).nnz_lower(),10*11/2)
    self.assertEqual(Sparsity.dense(10,10).nnz_upper(),10*11/2)
    self.assertEqual(Sparsity.dense(10,10).nnz_diag(),10)

    self.assertEqual(sparsify(DM([[1,1,0],[1,0,1],[0,0,0]])).nnz_diag(),1)
    self.assertEqual(sparsify(DM([[1,1,0],[1,0,1],[0,0,0]])).nnz_lower(),2)
    self.assertEqual(sparsify(DM([[1,1,0],[1,0,1],[0,0,0]])).nnz_upper(),3)

  def test_tril2symm(self):
    a = DM(Sparsity.upper(3),list(range(Sparsity.upper(3).nnz()))).T
    s = tril2symm(a)
    self.checkarray(s,DM([[0,1,3],[1,2,4],[3,4,5]]))

    with self.assertRaises(Exception):
      tril2symm(DM.ones(5,3))

    print(DM.ones(5,5).nnz_upper()-DM.ones(5,5).nnz_diag())

    with self.assertRaises(Exception):
      tril2symm(DM.ones(5,5))

  def test_not_null(self):
    x = MX.sym('x',3,1)
    sp = Sparsity.upper(2)
    MX(sp,x)

  def test_segfault(self):
    x = MX.sym('x',10,1)
    sp = Sparsity.upper(2)
    y = triu2symm(MX(sp,x[1:4]))
    f = Function("f", [x],[y])

  def test_append_empty(self):
    a = vertcat(DM(0,0),DM(0,2))

    self.assertEqual(a.size1(),0)
    self.assertEqual(a.size2(),2)

    a = vertcat(DM(0,0),DM(2,0),DM(3,0))

    self.assertEqual(a.size1(),5)
    self.assertEqual(a.size2(),0)

  def test_vertcat_empty(self):
    a = DM(0,2)
    v = vertcat(a,a)

    self.assertEqual(v.size1(),0)
    self.assertEqual(v.size2(),2)

    a = DM(2,0)
    v = vertcat(a,a)

    self.assertEqual(v.size1(),4)
    self.assertEqual(v.size2(),0)

  def test_vertsplit(self):
    a = DM(Sparsity.upper(5),list(range(int(5*6/2)))).T
    v = vertsplit(a,[0,2,4,5])

    self.assertEqual(len(v),3)
    self.checkarray(v[0],DM([[0,0,0,0,0],[1,2,0,0,0]]))
    self.checkarray(v[1],DM([[3,4,5,0,0],[6,7,8,9,0]]))
    self.checkarray(v[2],DM([[10,11,12,13,14]]))

    v = vertsplit(a)
    self.assertEqual(len(v),a.size1())
    self.checkarray(v[0],DM([[0,0,0,0,0]]))
    self.checkarray(v[1],DM([[1,2,0,0,0]]))
    self.checkarray(v[2],DM([[3,4,5,0,0]]))
    self.checkarray(v[3],DM([[6,7,8,9,0]]))
    self.checkarray(v[4],DM([[10,11,12,13,14]]))

    v = vertsplit(a,[0,2,4,5])
    self.assertEqual(len(v),3)
    self.checkarray(v[0],DM([[0,0,0,0,0],[1,2,0,0,0]]))
    self.checkarray(v[1],DM([[3,4,5,0,0],[6,7,8,9,0]]))
    self.checkarray(v[2],DM([[10,11,12,13,14]]))

    v = vertsplit(a,[0,0,3,5])
    self.assertEqual(len(v),3)
    self.assertEqual(v[0].size1(),0)
    self.assertEqual(v[0].size2(),5)
    self.checkarray(v[1],DM([[0,0,0,0,0],[1,2,0,0,0],[3,4,5,0,0]]))
    self.checkarray(v[2],DM([[6,7,8,9,0],[10,11,12,13,14]]))

  def test_horzsplit(self):
    a = DM(Sparsity.upper(5),list(range(int(5*6/2)))).T
    v = horzsplit(a,[0,2,4,5])

    self.assertEqual(len(v),3)
    self.checkarray(v[0],DM([[0,0],[1,2],[3,4],[6,7],[10,11]]))
    self.checkarray(v[1],DM([[0,0],[0,0],[5,0],[8,9],[12,13]]))
    self.checkarray(v[2],DM([[0],[0],[0],[0],[14]]))

    v = horzsplit(a)
    self.assertEqual(len(v),a.size1())
    self.checkarray(v[0],DM([0,1,3,6,10]))
    self.checkarray(v[1],DM([0,2,4,7,11]))
    self.checkarray(v[2],DM([0,0,5,8,12]))
    self.checkarray(v[3],DM([0,0,0,9,13]))
    self.checkarray(v[4],DM([0,0,0,0,14]))

    v = horzsplit(a,[0,2,4,5])
    self.assertEqual(len(v),3)
    self.checkarray(v[0],DM([[0,0],[1,2],[3,4],[6,7],[10,11]]))
    self.checkarray(v[1],DM([[0,0],[0,0],[5,0],[8,9],[12,13]]))
    self.checkarray(v[2],DM([[0],[0],[0],[0],[14]]))

    v = horzsplit(a,[0,0,3,5])
    self.assertEqual(len(v),3)
    self.assertEqual(v[0].size1(),5)
    self.assertEqual(v[0].size2(),0)
    self.checkarray(v[1],DM([[0,0,0],[1,2,0],[3,4,5],[6,7,8],[10,11,12]]))
    self.checkarray(v[2],DM([[0,0],[0,0],[0,0],[9,0],[13,14]]))

  def test_blocksplit(self):
    a = DM(Sparsity.upper(5),list(range(int(5*6/2)))).T
    v = blocksplit(a,[0,2,4,5],[0,1,3,5])

    self.checkarray(v[0][0],DM([0,1]))
    self.checkarray(v[0][1],DM([[0,0],[2,0]]))
    self.checkarray(v[1][0],DM([3,6]))
    self.checkarray(blockcat(v),a)

  def test_solve(self):
    import random

    spA = [
      Sparsity.dense(1,1)
    ]

    for n in range(2,5):
      spA+= [
        Sparsity.diag(n),
        Sparsity.dense(n,n),
        Sparsity.lower(n),
        Sparsity.lower(n).T,
        Sparsity.banded(n,1),
        diagcat(*[Sparsity.diag(n),Sparsity.dense(n,n)]),
        diagcat(*[Sparsity.diag(n),Sparsity.lower(n)]),
        diagcat(*[Sparsity.diag(n),Sparsity.lower(n).T]),
        diagcat(*[Sparsity.lower(n),Sparsity.lower(n).T]),
        Sparsity.diag(n)+Sparsity.rowcol([0],[n-1],n,n),
        Sparsity.diag(n)+Sparsity.rowcol([0,n-1],[n-1,0],n,n),
        Sparsity.diag(n)+Sparsity.triplet(n,n,[0],[n-1]),
        Sparsity.diag(n)+Sparsity.triplet(n,n,[0,n-1],[n-1,0]),
      ]

    for sA in spA:

      random.seed(1)
      a = DM(sA,[random.random() for i in range(sA.nnz())])
      A = SX.sym("a",a.sparsity())
      for sB in [ Sparsity.dense(a.size1(),1), vertcat(Sparsity.dense(1,1),Sparsity(a.size1()-1,1)),Sparsity.lower(a.size1()),Sparsity.lower(a.size1()).T]:

        b = DM(sB,[random.random() for i in range(sB.nnz())])
        B = SX.sym("B",b.sparsity())
        C = solve(A,B)

        f = Function("f", [A,B],[C])


        c = f(a,b)

        c_ref = DM(linalg.solve(a,b))
        c_ref = sparsify(c_ref)

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
    a = sparsify(DM([[1,0,6],[2,7,0]]))
    b = sparsify(DM([[1,0,0],[2,3,7],[0,0,9],[1,12,13]]))

    c_ = c.kron(a,b)

    self.assertEqual(c_.size1(),a.size1()*b.size1())
    self.assertEqual(c_.size2(),a.size2()*b.size2())
    self.assertEqual(c_.nnz(),a.nnz()*b.nnz())

    self.checkarray(c_,numpy.kron(a,b))

  def test_vec_kron(self):
    A = SX.sym("A",2,3)
    B = SX.sym("B",4,5)
    P = SX.sym("P",A.size2(),B.size1())

    f = Function("f", [vec(P.T),A,B],[vec(mtimes([A,P,B]).T)])

    J = f.jacobian_old(0, 0)
    J_in = []
    J_in.append(numpy.random.rand(*vec(P.T).shape))
    J_in.append(numpy.random.rand(*A.shape))
    J_in.append(numpy.random.rand(*B.shape))

    res, _  =  J(*J_in)

    ref =  kron(J_in[1],J_in[2].T)

    self.checkarray(res,ref)

  def test_repmat(self):
    a = DM([[1,2],[3,4],[5,6]])
    self.checkarray(repmat(a,2,3),kron(DM.ones(2,3),a))

  def test_triu(self):
    a = DM([[1,2],[3,4]])
    b = triu(a)
    self.checkarray(b, DM([[1,2],[0,4]]) )


  def test_tril(self):
    a = DM([[1,2],[3,4]])
    b = tril(a)
    self.checkarray(b, DM([[1,0],[3,4]]) )

  def test_nz(self):
    a = sparsify(IM([[1,2],[0,0],[3,4]]))
    self.checkarray(a.nz[:], IM([1,3,2,4]) )
    self.checkarray(len(a.nz), 4 )
    self.checkarray(a.nz[:-1], IM([1,3,2]) )
    self.checkarray(a.nz[0], IM([1]) )

  def test_norm_inf_mul(self):
    numpy.random.seed(0)

    A = numpy.random.random((10,2))
    B = numpy.random.random((2,8))

    self.checkarray(norm_inf_mul(A,B),norm_inf(mtimes(A,B)))
    self.checkarray(DM(norm_0_mul(A,B)),mtimes(A,B).nnz())

    # Sparse
    for i in range(5):
      A[numpy.random.randint(A.shape[0]),numpy.random.randint(A.shape[1])] = 0
      B[numpy.random.randint(B.shape[0]),numpy.random.randint(B.shape[1])] = 0

    A = sparsify(A)
    B = sparsify(B)

    self.checkarray(norm_inf_mul(A,B),norm_inf(mtimes(A,B)))
    self.checkarray(DM(norm_0_mul(A,B)),mtimes(A,B).nnz())


    A = numpy.random.random((8,2))
    B = numpy.random.random((2,10))

    self.checkarray(norm_inf_mul(A,B),norm_inf(mtimes(A,B)))
    self.checkarray(DM(norm_0_mul(A,B)),mtimes(A,B).nnz())

    # Sparse
    for i in range(5):
      A[numpy.random.randint(A.shape[0]),numpy.random.randint(A.shape[1])] = 0
      B[numpy.random.randint(B.shape[0]),numpy.random.randint(B.shape[1])] = 0

    A = sparsify(A)
    B = sparsify(B)

    self.checkarray(norm_inf_mul(A,B),norm_inf(mtimes(A,B)))
    self.checkarray(DM(norm_0_mul(A,B)),mtimes(A,B).nnz())

  def  test_mul3_issue_1465(self):
    with self.assertRaises(Exception):
      w = SX.sym("w",2,1)
      Q = np.eye(2)
      mtimes(w.T,Q,w)

  def test_chol(self):
    numpy.random.seed(0)

    for i in range(4):
      A = numpy.random.random((3,3))
      H = mtimes(A,A.T)

      R = chol(H)

      assert R.is_triu()
      self.checkarray(mtimes(R.T,R),H)
  def test_skew(self):
    x = DM([1,7,13])
    self.checkarray(inv_skew(skew(x)),x)
    y = DM([0.2,0.9,0.4])
    self.checkarray(mtimes(skew(x),y),cross(x,y))

  def test_nz_overflow(self):
    d = DM([2,3])
    r = d.nz[:]
    self.checkarray(r,d)

  def test_DMcrash(self):
    with self.assertRaises(Exception):
      DM([DM([1,2]),DM([1,2])])
    a = DM([DM([1]),DM([2])])
    self.checkarray(a,DM([1,2]))

  def test_sparsity_operation(self):
    L = [DM(1), DM(Sparsity(1,1),1), DM(Sparsity(2,1),1), DM(Sparsity.dense(2,1),1)]

    for a in L:
      for b in L:
        c = a*b

        if a.nnz()==0 or b.nnz()==0:
          self.assertTrue(c.nnz()==0)
        else:
          self.assertTrue(c.nnz()>0)

    self.assertTrue(sum2(IM(Sparsity(1,1),1)).nnz()==0)

  def test_matlab_operations(self):

    data = [ np.array([[1,3],[11,17]]) , np.array([[1,3]]) ,np.array([[1],[3]]), np.array([[3]])]

    for A in data:
      B = reshape(DM(A),A.shape)
      #self.checkarray(np.cumsum(A),cumsum(B))
      self.checkarray(np.cumsum(A,0),cumsum(B,0))
      self.checkarray(np.cumsum(A,1),cumsum(B,1))

      #self.checkarray(np.diff(A),diff(B))

      #self.checkarray(np.diff(A,1),diff(B,1))
      #self.checkarray(np.diff(A,2),diff(B,2))

      self.checkarray(np.diff(A,1,0),diff(B,1,0))
      #self.checkarray(np.diff(A,1,1),diff(B,1,1))


  def test_singular_repmat(self):
    for X in [DM, SX, MX, Sparsity]:
      for n_b in [0,2]:
        for m_b in [0,2]:
          b = X(n_b, m_b)

          for n in [0,3]:
            for m in [0,3]:
              self.assertEqual(repmat(b, n, m).shape,(n_b * n, m_b * m))

if __name__ == '__main__':
    unittest.main()
