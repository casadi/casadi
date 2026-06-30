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
import casadi as ca
import numpy as np
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

    a = ca.Sparsity(4,5)
    for i in nza:
      a.add_nz(i[0],i[1])

    b = ca.Sparsity(4,5)
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

    a = ca.Sparsity(4,5)
    for i in nza:
      a.add_nz(i[0],i[1])

    b = ca.Sparsity(4,5)
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

    a = ca.Sparsity(4,5)
    for i in nza:
      a.add_nz(i[0],i[1])

    A = ca.DM.ones(a)
    Ad = ca.DM(array(A))
    for i in a.find():
      self.assertEqual(Ad.nz[i],1)

  def test_find_nonzero(self):
    numpy.random.seed(0)
    d = self.randDM(20,10,0.6)
    sp = d.sparsity()

    sp2 = ca.Sparsity.nonzeros(20,10,sp.find())
    self.assertTrue(sp==sp2)

  def test_enlarge(self):
    self.message("enlarge")
    import numpy
    self.message(":dense")
    #sp = Sparsity(3,4,[1,2,1],[0,2,2,3])
    sp = ca.Sparsity.dense(3,4)

    col = [1,2,4]
    row = [0,3,4,6]
    sp.enlarge(7,8,col,row)

    z = numpy.zeros((7,8))
    for i in col:
      for j in row:
        z[i,j]=1

    self.checkarray(ca.DM.ones(sp),z,"enlarge")
    self.message(":sparse")
    sp = ca.Sparsity(4,3,[0,2,2,3],[1,2,1]).T
    n = ca.DM.ones(sp)
    z = numpy.zeros((7,8))
    for i in range(3):
      for j in range(4):
          z[col[i],row[j]]= n[i,j]
    sp.enlarge(7,8,[1,2,4],[0,3,4,6])

    self.checkarray(ca.DM.ones(sp),z,"enlarge")

  def tomatrix(self,s):
    d = ca.DM.ones(s)
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

    a = ca.Sparsity(4,5)
    for i in nza:
      a.add_nz(i[0],i[1])

    b = ca.Sparsity.triplet(4,5,[i[0] for i in nza],[i[1] for i in nza])
    self.checkarray(self.tomatrix(a),self.tomatrix(b),"rowcol")

  def test_rowcol(self):
    self.message("rowcol constructor")

    r = [0,1,3]
    c = [1,4]
    a = ca.Sparsity(4,5)
    for i in r:
      for j in c:
        a.add_nz(i,j)

    b = ca.Sparsity.rowcol(r,c,4,5)
    self.checkarray(self.tomatrix(a),self.tomatrix(b),"rowcol")

  def test_reshape(self):
    self.message("Reshape")
    nza = set([  (0,0),
             (0,1),
             (2,0),
             (2,3),
             (2,4),
             (3,1)])

    a = ca.Sparsity(4,5)
    for i in nza:
      a.add_nz(i[0],i[1])

    A=self.tomatrix(a).full()
    B=self.tomatrix(ca.reshape(a,2,10)).full()
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

    a = ca.Sparsity(4,5)
    for i in nza:
      a.add_nz(i[0],i[1])

    A=self.tomatrix(a).full()
    B=self.tomatrix(ca.vec(a)).full()
    B_=numpy.reshape(A,(20,1))

    self.checkarray(B,B_,"reshape")


  def test_refcount(self):
      x = ca.DM(ca.Sparsity.lower(4),5)
      s = (x @ x).sparsity()
      self.assertEqual(s.numel(),16)

  def test_splower(self):
    sp = ca.Sparsity(4,3,[0,2,2,3],[1,2,1])
    print(array(sp))
    print(array(ca.tril(sp)))
    print(sp.get_lower())


  def test_diag(self):
    self.message("diag")
    A = ca.Sparsity(5,5)
    A.add_nz(1,1)
    A.add_nz(2,4)
    A.add_nz(3,3)

    sp, mapping = A.get_diag()
    B = ca.DM.ones(sp)

    self.checkarray(array([[0],[1],[0],[1],[0]]),B,"get_diag(matrix)")
    self.checkarray(array([0,1]),array(list(mapping)),"get_diag(vector)")

    #print B

    A = ca.Sparsity(5,1)
    A.add_nz(1,0)
    A.add_nz(2,0)
    A.add_nz(4,0)

    sp, mapping = A.get_diag()
    B = ca.DM.ones(sp)

    self.checkarray(array([[0,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,0,0],[0,0,0,0,1]]),B,"get_diag(vector)")

    self.checkarray(array([0,1,2]),array(list(mapping)),"get_diag(vector)")

    A = ca.Sparsity(1,5)
    A.add_nz(0,1)
    A.add_nz(0,2)
    A.add_nz(0,4)

    sp, mapping = A.get_diag()
    B = ca.DM.ones(sp)

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

    a = ca.Sparsity(5,5)
    for i in nza:
      a.add_nz(i[0],i[1])

    b = ca.MX.sym("b",a)

    self.assertRaises(Exception,lambda: b[ca.Sparsity.diag(3)])

    d = ca.Sparsity.diag(5)
    c = b[d]

    self.assertTrue(c.sparsity()==d)

    f = ca.Function('f', [b],[c])
    fin = ca.DM(b.sparsity(),list(range(1,len(nza)+1)))
    f_out = f(fin)

    self.checkarray(ca.DM(f_out.nonzeros()),ca.DM([1,0,0,7,0]),"sparsity index")

  def test_get_ccs(self):
    self.message("CCS format")
    nza = set([  (0,0),
             (0,1),
             (2,0),
             (2,3),
             (3,3),
             (2,4),
             (3,1)])

    a = ca.Sparsity(4,5)
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
      A = ca.diagcat(*Ai)

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
      A = ca.diagcat(*Ai)

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

    A = ca.DM(6,4)
    A[0,0] = 1
    A[1,2] = 1
    A[2,2] = 1
    A[5,3] = 1

    ret, rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock = A.sparsity().btf()

    # Checked with CSparse
    self.checkarray(ca.DM([ret]),ca.DM([4]))
    self.checkarray(rowperm,ca.DM([2, 3, 4, 1, 0, 5]).T)
    self.checkarray(colperm,ca.DM([ 2,0,3,1]).T)
    self.checkarray(rowblock,ca.DM([ 0, 4,5,6,6]).T)
    self.checkarray(colblock,ca.DM([ 0, 1,2,3,4]).T)
    self.checkarray(coarse_rowblock,ca.DM([ 0, 3,4,6,6]).T)
    self.checkarray(coarse_colblock,ca.DM([ 0, 1,3,3,4]).T)


    A = ca.DM(6,4)
    A[0,0] = 1
    A[1,2] = 1
    A[2,2] = 1
    A[5,3] = 1
    A[4,1] = 1
    A[3,0] = 1

    A.sparsity().spy()

    ret, rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock = A.sparsity().btf()

    # Checked with CSparse
    self.checkarray(ca.DM([ret]),ca.DM([3]))
    self.checkarray(rowperm,ca.DM([2,3,0,1,4,5]).T)
    self.checkarray(colperm,ca.DM([ 0, 2, 1, 3]).T)
    self.checkarray(rowblock,ca.DM([ 0, 4,5,6]).T)
    self.checkarray(colblock,ca.DM([ 0, 2,3,4]).T)
    self.checkarray(coarse_rowblock,ca.DM([ 0, 2, 4,6,6]).T)
    self.checkarray(coarse_colblock,ca.DM([ 0, 2,4,4,4]).T)

    A = ca.DM(6,4)
    A[0,0] = 1
    A[1,2] = 1
    A[2,2] = 1
    A[5,3] = 1
    A[4,1] = 1
    A[3,0] = 1
    A = A + ca.DM.eye(6)[:,:4]

    A.sparsity().spy()

    ret, rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock = A.sparsity().btf()

    # Checked with CSparse
    self.checkarray(ca.DM([ret]),ca.DM([1]))
    self.checkarray(rowperm,ca.DM([4, 5, 0, 1, 2, 3]).T)
    self.checkarray(colperm,ca.DM([ 0, 1, 2, 3]).T)
    self.checkarray(rowblock,ca.DM([ 0, 6]).T)
    self.checkarray(colblock,ca.DM([ 0, 4]).T)
    self.checkarray(coarse_rowblock,ca.DM([ 0, 2, 6,6,6]).T)
    self.checkarray(coarse_colblock,ca.DM([ 0, 4,4,4,4]).T)


  def test_jacsparsityHierarchical(self):

    X = ca.SX.sym("X",100)
    P = ca.SX.sym("P",1000)

    optvar = ca.vertcat(*[X,P])

    p = ca.SX.sym("p")

    g = ca.Function('g', [optvar,p],[X*p], {'verbose':True})

    J = jacobian_old(g, 0, 0)

    self.assertTrue(ca.DM(J.sparsity_out(0))[:,:X.nnz()].sparsity()==ca.Sparsity.diag(100))

    X = ca.SX.sym("X",100)
    P = ca.SX.sym("P",1000)

    p = ca.SX.sym("p")

    g = ca.Function('g', [X,P,p],[ca.vertcat(*[X*p,P])], {'verbose':True})

    J = jacobian_old(g, 0, 0)

    self.assertTrue(ca.DM(J.sparsity_out(0))[:X.nnz(),:].sparsity()==ca.Sparsity.diag(100))

  @memory_heavy()
  def test_jacsparsityHierarchicalSymm(self):
    ca.GlobalOptions.setHierarchicalSparsity(False)
    sp = ca.Sparsity.banded(4129,1)

    x = ca.MX.sym("x",sp.size1())

    H = ca.bilin(ca.MX.ones(sp),x,x)
    sp2 = ca.hessian(H,x)[0].sparsity()
    self.assertTrue(sp==sp2)


  def test_rowcol(self):
    n = 3

    s = ca.Sparsity.rowcol([n-1,0],[0,n-1],n,n)
    self.checkarray(s.colind(),[0,2,2,4])
    self.checkarray(s.row(),[0,2,0,2])

  def test_inverse(self):
    numpy.random.seed(0)
    d = self.randDM(20,20,0.6,symm=True)
    sp = d.sparsity()

    for sp in [sp,ca.Sparsity.dense(4,4),ca.Sparsity(4,4),ca.Sparsity.lower(4),ca.Sparsity.lower(4).T]:

      d = ca.DM.ones(sp)

      dt = ca.sparsify(1-d)
      dt = ca.DM.ones(dt.sparsity())

      trial = ca.DM.ones(sp.pattern_inverse())

      d.print_dense()
      dt.print_dense()
      trial.print_dense()

      self.checkarray(trial,dt)

  def test_kron(self):
    a = ca.sparsify(ca.DM([[1,0,6],[2,7,0]]))
    b = ca.sparsify(ca.DM([[1,0,0],[2,3,7],[0,0,9],[1,12,13]]))

    c_ = c.kron(a.sparsity(),b.sparsity())

    self.assertEqual(c_.size1(),a.size1()*b.size1())
    self.assertEqual(c_.size2(),a.size2()*b.size2())
    self.assertEqual(c_.nnz(),a.nnz()*b.nnz())

    self.checkarray(ca.DM(c_,1),ca.DM(c.kron(a,b).sparsity(),1))

  def test_nz_method(self):
    n = 20
    m = 25
    import random
    random.seed(0)
    numpy.random.seed(0)
    d = self.randDM(n,m,0.5)
    D = ca.densify(ca.vec(d))
    dn = ca.DM(d.nonzeros())
    sp = d.sparsity()
    z = np.unique([random.randint(0,n*m-1) for i in range(200)])
    zres = sp.get_nz(z)
    A = dn[[e for e in zres if e>=0]]
    B = D[[e for e,k in zip(z,zres) if k>=0]]
    self.checkarray(A,B)
    self.assertFalse(np.any(D[[e for e,k in zip(z,zres) if k==-1]]))

  def test_serialize(self):
    for a in [ca.Sparsity(), ca.Sparsity.dense(4,5), ca.Sparsity.lower(5)]:
      b = ca.Sparsity.deserialize(a.serialize())
      if a.is_null():
        self.assertTrue(b.is_null())
      else:
        self.checkarray(ca.DM(a,1),ca.DM(b,1))

  def test_is_subset(self):

      pairs = [ (ca.Sparsity.lower(3), ca.Sparsity.dense(3,3)),
                (ca.Sparsity.diag(3), ca.Sparsity.dense(3,3)),
                (ca.Sparsity.diag(3), ca.Sparsity.lower(3)),
                (ca.Sparsity(3,3), ca.Sparsity.lower(3)),
      ]

      for L,R in pairs:
        self.assertTrue(L.is_subset(R))
        self.assertFalse(R.is_subset(L))

  def test_is_compactible(self):
    self.message("is_compactible")
    # A sparsity is compactible iff its nonzero pattern is the Cartesian
    # product of a row subset and a column subset; the CCS buffer then
    # equals a column-major dense row.size() x col.size() matrix.

    # Cartesian (0,1), (2,1), (0,3), (2,3) -- compactible with row {0,2}, col {1,3}
    sp = ca.Sparsity(4, 5)
    sp.add_nz(0, 1); sp.add_nz(2, 1)
    sp.add_nz(0, 3); sp.add_nz(2, 3)
    is_compact, row_support, col_support = sp.is_compactible()
    self.assertTrue(is_compact)
    self.checkarray(row_support, [0, 2])
    self.checkarray(col_support, [1, 3])

    # Drop one nonzero: not a Cartesian product anymore.
    sp_partial = ca.Sparsity(4, 5)
    sp_partial.add_nz(0, 1); sp_partial.add_nz(2, 1); sp_partial.add_nz(2, 3)
    self.assertFalse(sp_partial.is_compactible()[0])

    # rowcol() builds a compactible pattern by construction.
    sp_rc = ca.Sparsity.rowcol([0, 2], [1, 3], 4, 5)
    is_compact, rs, cs = sp_rc.is_compactible()
    self.assertTrue(is_compact)
    self.checkarray(rs, [0, 2])
    self.checkarray(cs, [1, 3])

    # Empty sparsity: trivially compactible with empty row/col support.
    sp_empty = ca.Sparsity(3, 4)
    is_compact, rs, cs = sp_empty.is_compactible()
    self.assertTrue(is_compact)
    self.checkarray(rs, [])
    self.checkarray(cs, [])

    # Dense sparsity: row support = all rows, col support = all cols.
    sp_dense = ca.Sparsity.dense(3, 4)
    is_compact, rs, cs = sp_dense.is_compactible()
    self.assertTrue(is_compact)
    self.checkarray(rs, [0, 1, 2])
    self.checkarray(cs, [0, 1, 2, 3])

    # Diagonal of a 4x4: nonempty rows/cols both = {0..3} but nnz=4 != 16.
    sp_diag = ca.Sparsity.diag(4)
    self.assertFalse(sp_diag.is_compactible()[0])


if __name__ == '__main__':
    unittest.main()
