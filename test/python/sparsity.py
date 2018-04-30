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
import unittest
from types import *
from helpers import *
import numpy
import random

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

    a = Sparsity(4,5)
    for i in nza:
      a.add_nz(i[0],i[1])

    b = Sparsity(4,5)
    for i in nzb:
      b.add_nz(i[0],i[1])

    c =a.unite(b)

    c = a + b
    self.assertEqual(c.nnz(),len(nza.union(nzb)))
    for k in range(c.nnz()):
      ind = (c.row(k),c.get_col()[k])
      self.assertTrue(ind in nza or ind in nzb)

  def test_intersection(self):
    self.message("Sparsity intersection")
    nza = set([  (0,0),
             (0,1),
             (2,0),
             (3,1),
             (2,3)])
    nzb = set([  (0,2),
             (0,0),
             (2,2),
             (2,3)])

    a = Sparsity(4,5)
    for i in nza:
      a.add_nz(i[0],i[1])

    b = Sparsity(4,5)
    for i in nzb:
      b.add_nz(i[0],i[1])

    c=a.intersect(b)
    for k in range(c.nnz()):
      ind = (c.row(k),c.get_col()[k])
      self.assertTrue(ind in nza and ind in nzb)

    c = a * b
    self.assertEqual(c.nnz(),len(nza.intersection(nzb)))
    for k in range(c.nnz()):
      ind = (c.row(k),c.get_col()[k])
      self.assertTrue(ind in nza and ind in nzb)

  def test_get_nzDense(self):
    self.message("get_nzDense")
    nza = set([  (0,0),(0,1),(2,0),(3,1)])

    a = Sparsity(4,5)
    for i in nza:
      a.add_nz(i[0],i[1])

    A = DM.ones(a)
    Ad = DM(array(A))
    for i in a.find():
      self.assertEqual(Ad.nz[i],1)
      
  def test_find_nonzero(self):
    numpy.random.seed(0)
    d = self.randDM(20,10,0.6)
    sp = d.sparsity()
    
    sp2 = Sparsity.nonzeros(20,10,sp.find())
    self.assertTrue(sp==sp2)

  def test_enlarge(self):
    self.message("enlarge")
    import numpy
    self.message(":dense")
    #sp = Sparsity(3,4,[1,2,1],[0,2,2,3])
    sp = Sparsity.dense(3,4)

    col = [1,2,4]
    row = [0,3,4,6]
    sp.enlarge(7,8,col,row)

    z = numpy.zeros((7,8))
    for i in col:
      for j in row:
        z[i,j]=1

    self.checkarray(DM.ones(sp),z,"enlarge")
    self.message(":sparse")
    sp = Sparsity(4,3,[0,2,2,3],[1,2,1]).T
    n = DM.ones(sp)
    z = numpy.zeros((7,8))
    for i in range(3):
      for j in range(4):
          z[col[i],row[j]]= n[i,j]
    sp.enlarge(7,8,[1,2,4],[0,3,4,6])

    self.checkarray(DM.ones(sp),z,"enlarge")

  def tomatrix(self,s):
    d = DM.ones(s)
    for k in range(d.nnz()):
      d.nz[k] = k+1
    return d

  def test_NZ(self):
    self.message("NZ constructor")
    nza = [  (0,0),
             (0,1),
             (2,0),
             (2,3),
             (2,4),
             (3,1)]

    a = Sparsity(4,5)
    for i in nza:
      a.add_nz(i[0],i[1])

    b = Sparsity.triplet(4,5,[i[0] for i in nza],[i[1] for i in nza])
    self.checkarray(self.tomatrix(a),self.tomatrix(b),"rowcol")

  def test_rowcol(self):
    self.message("rowcol constructor")

    r = [0,1,3]
    c = [1,4]
    a = Sparsity(4,5)
    for i in r:
      for j in c:
        a.add_nz(i,j)

    b = Sparsity.rowcol(r,c,4,5)
    self.checkarray(self.tomatrix(a),self.tomatrix(b),"rowcol")

  def test_reshape(self):
    self.message("Reshape")
    nza = set([  (0,0),
             (0,1),
             (2,0),
             (2,3),
             (2,4),
             (3,1)])

    a = Sparsity(4,5)
    for i in nza:
      a.add_nz(i[0],i[1])

    A=self.tomatrix(a).full()
    B=self.tomatrix(casadi.reshape(a,2,10)).full()
    B_=numpy.reshape(A.T,(10,2)).T

    self.checkarray(B,B_,"reshape")

  def test_vec(self):
    return # This test doesn't make much sense
    self.message("vec")
    nza = set([  (0,0),
             (0,1),
             (2,0),
             (2,3),
             (2,4),
             (3,1)])

    a = Sparsity(4,5)
    for i in nza:
      a.add_nz(i[0],i[1])

    A=self.tomatrix(a).full()
    B=self.tomatrix(vec(a)).full()
    B_=numpy.reshape(A,(20,1))

    self.checkarray(B,B_,"reshape")


  def test_refcount(self):
      x = DM(Sparsity.lower(4),5)
      s = mtimes(x,x).sparsity()
      self.assertEqual(s.numel(),16)

  def test_splower(self):
    sp = Sparsity(4,3,[0,2,2,3],[1,2,1])
    print(array(sp))
    print(array(tril(sp)))
    print(sp.get_lower())


  def test_diag(self):
    self.message("diag")
    A = Sparsity(5,5)
    A.add_nz(1,1)
    A.add_nz(2,4)
    A.add_nz(3,3)

    sp, mapping = A.get_diag()
    B = DM.ones(sp)

    self.checkarray(array([[0],[1],[0],[1],[0]]),B,"get_diag(matrix)")
    self.checkarray(array([0,1]),array(list(mapping)),"get_diag(vector)")

    #print B

    A = Sparsity(5,1)
    A.add_nz(1,0)
    A.add_nz(2,0)
    A.add_nz(4,0)

    sp, mapping = A.get_diag()
    B = DM.ones(sp)

    self.checkarray(array([[0,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,0,0],[0,0,0,0,1]]),B,"get_diag(vector)")

    self.checkarray(array([0,1,2]),array(list(mapping)),"get_diag(vector)")

    A = Sparsity(1,5)
    A.add_nz(0,1)
    A.add_nz(0,2)
    A.add_nz(0,4)

    sp, mapping = A.get_diag()
    B = DM.ones(sp)

    self.checkarray(array([[0,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,0,0],[0,0,0,0,1]]),B,"get_diag(vector)")

    self.checkarray(array([0,1,2]),array(list(mapping)),"get_diag(vector)")

  def test_sparsityindex(self):
    self.message("sparsity indexing")
    nza = set([  (0,0),
             (0,1),
             (2,0),
             (2,3),
             (3,3),
             (2,4),
             (3,1),
             (4,1)])

    a = Sparsity(5,5)
    for i in nza:
      a.add_nz(i[0],i[1])

    b = MX.sym("b",a)

    self.assertRaises(Exception,lambda: b[Sparsity.diag(3)])

    d = Sparsity.diag(5)
    c = b[d]

    self.assertTrue(c.sparsity()==d)

    f = Function('f', [b],[c])
    fin = DM(b.sparsity(),list(range(1,len(nza)+1)))
    f_out = f(fin)

    self.checkarray(DM(f_out.nonzeros()),DM([1,0,0,7,0]),"sparsity index")

  def test_get_ccs(self):
    self.message("CCS format")
    nza = set([  (0,0),
             (0,1),
             (2,0),
             (2,3),
             (3,3),
             (2,4),
             (3,1)])

    a = Sparsity(4,5)
    for i in nza:
      a.add_nz(i[0],i[1])

    A1, B1= a.get_ccs()

    A2, B2 = (a.T).get_crs()

    print(A1, B1)
    print(A2, B2)

  def test_dm_diagcat_dense(self):
    self.message("Dulmage-Mendelsohn")
    random.seed(0)
    numpy.random.seed(0)
    for k in range(20):
      Ai = [self.randDM(d,d,1) for i,d in enumerate ([random.randint(1,10) for j in range(10)])]
      A = diagcat(*Ai)

      #A.sparsity().spy()
      perm =  numpy.random.permutation(list(range(A.size1())))

      AP = A[perm,perm]
      #AP.sparsity().spy()

      ret, rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock = AP.sparsity().btf()

      Ar = AP[rowperm,colperm]

      ST = Ar.sparsity()

      blocks = []
      acc = -1
      mc = 0
      for i in range(0,Ar.size1()):
        mc = max(ST.row()[ST.colind()[i+1]-1],mc)
        if mc==i:
          blocks.append(i-acc)
          acc = i

      truth = [i.size1() for i in Ai]
      tryme = blocks

      truth.sort()
      tryme.sort()

      self.checkarray(truth,tryme)

  def test_scc_diagcat_sparse(self):
    self.message("scc")
    random.seed(0)
    numpy.random.seed(0)
    for k in range(20):
      Ai = [self.randDM(d,d,0.6,symm=True) for i,d in enumerate ([random.randint(1,10) for j in range(10)])]
      A = diagcat(*Ai)

      #A.sparsity().spy()
      perm =  numpy.random.permutation(list(range(A.size1())))

      AP = A[perm,perm]
      #AP.sparsity().spy()

      n,p,r = AP.sparsity().scc()

      Ar = AP[p,p]

      #print "permute"
      #Ar.sparsity().spy()

      ST = Ar.sparsity()

      blocks = []
      acc = -1
      mc = 0
      for i in range(0,Ar.size1()):
        mc = max(ST.row()[ST.colind()[i+1]-1],mc)
        if mc==i:
          blocks.append(i-acc)
          acc = i

      truth = [i.size1() for i in Ai]
      tryme = blocks

      self.assertTrue(n>=len(truth))
      self.assertTrue(n>=len(tryme))

  def test_dm(self):

    A = DM(6,4)
    A[0,0] = 1
    A[1,2] = 1
    A[2,2] = 1
    A[5,3] = 1

    ret, rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock = A.sparsity().btf()

    # Checked with CSparse
    self.checkarray(DM([ret]),DM([4]))
    self.checkarray(rowperm,DM([2, 3, 4, 1, 0, 5]).T)
    self.checkarray(colperm,DM([ 2,0,3,1]).T)
    self.checkarray(rowblock,DM([ 0, 4,5,6,6]).T)
    self.checkarray(colblock,DM([ 0, 1,2,3,4]).T)
    self.checkarray(coarse_rowblock,DM([ 0, 3,4,6,6]).T)
    self.checkarray(coarse_colblock,DM([ 0, 1,3,3,4]).T)


    A = DM(6,4)
    A[0,0] = 1
    A[1,2] = 1
    A[2,2] = 1
    A[5,3] = 1
    A[4,1] = 1
    A[3,0] = 1

    A.sparsity().spy()

    ret, rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock = A.sparsity().btf()

    # Checked with CSparse
    self.checkarray(DM([ret]),DM([3]))
    self.checkarray(rowperm,DM([2,3,0,1,4,5]).T)
    self.checkarray(colperm,DM([ 0, 2, 1, 3]).T)
    self.checkarray(rowblock,DM([ 0, 4,5,6]).T)
    self.checkarray(colblock,DM([ 0, 2,3,4]).T)
    self.checkarray(coarse_rowblock,DM([ 0, 2, 4,6,6]).T)
    self.checkarray(coarse_colblock,DM([ 0, 2,4,4,4]).T)

    A = DM(6,4)
    A[0,0] = 1
    A[1,2] = 1
    A[2,2] = 1
    A[5,3] = 1
    A[4,1] = 1
    A[3,0] = 1
    A = A + DM.eye(6)[:,:4]

    A.sparsity().spy()

    ret, rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock = A.sparsity().btf()

    # Checked with CSparse
    self.checkarray(DM([ret]),DM([1]))
    self.checkarray(rowperm,DM([4, 5, 0, 1, 2, 3]).T)
    self.checkarray(colperm,DM([ 0, 1, 2, 3]).T)
    self.checkarray(rowblock,DM([ 0, 6]).T)
    self.checkarray(colblock,DM([ 0, 4]).T)
    self.checkarray(coarse_rowblock,DM([ 0, 2, 6,6,6]).T)
    self.checkarray(coarse_colblock,DM([ 0, 4,4,4,4]).T)


  def test_jacsparsityHierarchical(self):

    X = SX.sym("X",100)
    P = SX.sym("P",1000)

    optvar = vertcat(*[X,P])

    p = SX.sym("p")

    g = Function('g', [optvar,p],[X*p], {'verbose':True})

    J = g.jacobian_old(0, 0)

    self.assertTrue(DM(J.sparsity_out(0))[:,:X.nnz()].sparsity()==Sparsity.diag(100))

    X = SX.sym("X",100)
    P = SX.sym("P",1000)

    p = SX.sym("p")

    g = Function('g', [X,p],[vertcat(*[X*p,P])], {'verbose':True})

    J = g.jacobian_old(0, 0)

    self.assertTrue(DM(J.sparsity_out(0))[:X.nnz(),:].sparsity()==Sparsity.diag(100))

  @memory_heavy()
  def test_jacsparsityHierarchicalSymm(self):
    GlobalOptions.setHierarchicalSparsity(False)
    sp = Sparsity.banded(4129,1)

    x = MX.sym("x",sp.size1())

    H = bilin(MX.ones(sp),x,x)
    sp2 = hessian(H,x)[0].sparsity()
    self.assertTrue(sp==sp2)


  def test_rowcol(self):
    n = 3

    s = Sparsity.rowcol([n-1,0],[0,n-1],n,n)
    self.checkarray(IM(s.colind()),IM([0,2,2,4]))
    self.checkarray(IM(s.row()),IM([0,2,0,2]))

  def test_inverse(self):
    numpy.random.seed(0)
    d = self.randDM(20,20,0.6,symm=True)
    sp = d.sparsity()

    for sp in [sp,Sparsity.dense(4,4),Sparsity(4,4),Sparsity.lower(4),Sparsity.lower(4).T]:

      d = IM.ones(sp)

      dt = sparsify(1-d)
      dt = IM.ones(dt.sparsity())

      trial = IM.ones(sp.pattern_inverse())

      d.print_dense()
      dt.print_dense()
      trial.print_dense()

      self.checkarray(trial,dt)

  def test_kron(self):
    a = sparsify(DM([[1,0,6],[2,7,0]]))
    b = sparsify(DM([[1,0,0],[2,3,7],[0,0,9],[1,12,13]]))

    c_ = c.kron(a.sparsity(),b.sparsity())

    self.assertEqual(c_.size1(),a.size1()*b.size1())
    self.assertEqual(c_.size2(),a.size2()*b.size2())
    self.assertEqual(c_.nnz(),a.nnz()*b.nnz())

    self.checkarray(IM(c_,1),IM(c.kron(a,b).sparsity(),1))

  def test_nz_method(self):
    n = 20
    m = 25
    import random
    random.seed(0)
    numpy.random.seed(0)
    d = self.randDM(n,m,0.5)
    D = densify(vec(d))
    dn = DM(d.nonzeros())
    sp = d.sparsity()
    z = np.unique([random.randint(0,n*m-1) for i in range(200)])
    zres = sp.get_nz(z)
    A = dn[[e for e in zres if e>=0]]
    B = D[[e for e,k in zip(z,zres) if k>=0]]
    self.checkarray(A,B)
    self.assertFalse(np.any(D[[e for e,k in zip(z,zres) if k==-1]]))    

if __name__ == '__main__':
    unittest.main()
