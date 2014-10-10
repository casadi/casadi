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
from numpy import *
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
    
    a = Sparsity.sparse(4,5)
    for i in nza:
      a.getNZ(i[0],i[1])
      
    b = Sparsity.sparse(4,5)  
    for i in nzb:
      b.getNZ(i[0],i[1])
      
    c,w =a.patternUnion(b)
    self.assertEquals(len(w),len(nza.union(nzb)))
    for k in range(len(w)):
      ind = (c.row(k),c.getCol()[k])
      if (ind in nza and ind in nzb):
        self.assertEquals(w[k],1 | 2)
      elif (ind in nza):
        self.assertEquals(w[k],1)
      elif (ind in nzb):
        self.assertEquals(w[k],2)
        
    c = a + b
    self.assertEquals(c.size(),len(nza.union(nzb)))
    for k in range(c.size()):
      ind = (c.row(k),c.getCol()[k])
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
    
    a = Sparsity.sparse(4,5)
    for i in nza:
      a.getNZ(i[0],i[1])
      
    b = Sparsity.sparse(4,5)  
    for i in nzb:
      b.getNZ(i[0],i[1])
    
    c,_=a.patternIntersection(b)
    for k in range(c.size()):
      ind = (c.row(k),c.getCol()[k])
      self.assertTrue(ind in nza and ind in nzb)
        
    c = a * b
    self.assertEquals(c.size(),len(nza.intersection(nzb)))
    for k in range(c.size()):
      ind = (c.row(k),c.getCol()[k])
      self.assertTrue(ind in nza and ind in nzb)
       
  def test_getNZDense(self):
    self.message("getNZDense")
    nza = set([  (0,0),(0,1),(2,0),(3,1)])
    
    a = Sparsity.sparse(4,5)
    for i in nza:
      a.getNZ(i[0],i[1])
      
    A = DMatrix(a,1)
    Ad = DMatrix(array(A))
    for i in a.getElements():
      self.assertEqual(Ad.nz[i],1)

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

    self.checkarray(DMatrix(sp,1),z,"enlarge")
    self.message(":sparse")
    sp = Sparsity(4,3,[0,2,2,3],[1,2,1]).T
    n = DMatrix(sp,1)
    z = numpy.zeros((7,8))
    for i in range(3):
      for j in range(4):
          z[col[i],row[j]]= n[i,j]
    sp.enlarge(7,8,[1,2,4],[0,3,4,6])
    
    self.checkarray(DMatrix(sp,1),z,"enlarge")
    
  def tomatrix(self,s):
    d = DMatrix(s,1)
    for k in range(d.size()):
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
    
    a = Sparsity.sparse(4,5)
    for i in nza:
      a.getNZ(i[0],i[1])
      
    b = Sparsity.triplet(4,5,[i[0] for i in nza],[i[1] for i in nza])
    self.checkarray(self.tomatrix(a),self.tomatrix(b),"rowcol")

  def test_rowcol(self):
    self.message("rowcol constructor")
    
    r = [0,1,3]
    c = [1,4]
    a = Sparsity.sparse(4,5)
    for i in r:
      for j in c:
        a.getNZ(i,j)
      
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
    
    a = Sparsity.sparse(4,5)
    for i in nza:
      a.getNZ(i[0],i[1])
      
    A=self.tomatrix(a).toArray()
    B=self.tomatrix(casadi.reshape(a,2,10)).toArray()
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
    
    a = Sparsity.sparse(4,5)
    for i in nza:
      a.getNZ(i[0],i[1])
      
    A=self.tomatrix(a).toArray()
    B=self.tomatrix(vec(a)).toArray()
    B_=numpy.reshape(A,(20,1))
    
    self.checkarray(B,B_,"reshape")
    
    
  def test_refcount(self):
      x = DMatrix(Sparsity.tril(4),5)
      s = mul(x,x).sparsity()
      self.assertEqual(s.numel(),16)
      
  def test_splower(self):
    sp = Sparsity(4,3,[0,2,2,3],[1,2,1])
    print array(sp)
    print array(tril(sp))
    print sp.getLowerNZ()
    
    
  def test_diag(self):
    self.message("diag")
    A = Sparsity.sparse(5,5)
    A.getNZ(1,1)
    A.getNZ(2,4)
    A.getNZ(3,3)
    
    sp, mapping = A.getDiag()
    B = DMatrix(sp,1)
    
    self.checkarray(array([[0],[1],[0],[1],[0]]),B,"getDiag(matrix)")
    self.checkarray(array([0,1]),array(list(mapping)),"getDiag(vector)")
    
    #print B
    
    A = Sparsity.sparse(5,1)
    A.getNZ(1,0)
    A.getNZ(2,0)
    A.getNZ(4,0)
    
    sp, mapping = A.getDiag()
    B = DMatrix(sp,1)
    
    self.checkarray(array([[0,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,0,0],[0,0,0,0,1]]),B,"getDiag(vector)")
    
    self.checkarray(array([0,1,2]),array(list(mapping)),"getDiag(vector)")
    
    A = Sparsity.sparse(1,5)
    A.getNZ(0,1)
    A.getNZ(0,2)
    A.getNZ(0,4)
    
    sp, mapping = A.getDiag()
    B = DMatrix(sp,1)
    
    self.checkarray(array([[0,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,0,0],[0,0,0,0,1]]),B,"getDiag(vector)")
    
    self.checkarray(array([0,1,2]),array(list(mapping)),"getDiag(vector)")
            
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
    
    a = Sparsity.sparse(5,5)
    for i in nza:
      a.getNZ(i[0],i[1])
      
    b = SX.sym("b",a)
    
    self.assertRaises(Exception,lambda: b[Sparsity.diag(3)])
    
    d = Sparsity.diag(5)
    c = b[d]

    self.assertTrue(c.sparsity()==d)
    
    f = SXFunction([b],[c])
    f.init()
    f.setInput(range(1,len(nza)+1))
    f.evaluate()
    
    self.checkarray(DMatrix(f.getOutput().data()),DMatrix([1,0,0,7,0]),"sparsity index")
    
    self.assertTrue(f.getOutput().data()[1]==0)
    
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
    
    a = Sparsity.sparse(5,5)
    for i in nza:
      a.getNZ(i[0],i[1])
      
    b = MX.sym("b",a)
    
    self.assertRaises(Exception,lambda: b[Sparsity.diag(3)])
    
    d = Sparsity.diag(5)
    c = b[d]

    self.assertTrue(c.sparsity()==d)
    
    f = MXFunction([b],[c])
    f.init()
    f.setInput(range(1,len(nza)+1))
    f.evaluate()
    
    self.checkarray(DMatrix(f.getOutput().data()),DMatrix([1,0,0,7,0]),"sparsity index")
    
  def test_getCCS(self):
    self.message("CCS format")
    nza = set([  (0,0),
             (0,1),
             (2,0),
             (2,3),
             (3,3),
             (2,4),
             (3,1)])
    
    a = Sparsity.sparse(4,5)
    for i in nza:
      a.getNZ(i[0],i[1])
      
    A1, B1= a.getCCS()
    
    A2, B2 = (a.T).getCRS()
    
    print A1, B1
    print A2, B2
    
  def test_dm_blkdiag_dense(self):
    self.message("Dulmage-Mendelsohn")
    random.seed(0)
    numpy.random.seed(0)
    for k in range(20):
      Ai = [self.randDMatrix(d,d,1) for i,d in enumerate ([random.randint(1,10) for j in range(10)])]
      A = blkdiag(Ai)
      
      #A.sparsity().spy()
      perm =  numpy.random.permutation(range(A.size1()))
      
      AP = A[perm,perm]
      #AP.sparsity().spy()

      ret, rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock = AP.sparsity().dulmageMendelsohn()

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

  def test_scc_blkdiag_sparse(self):
    self.message("stronglyConnectedComponents")
    random.seed(0)
    numpy.random.seed(0)
    for k in range(20):
      Ai = [self.randDMatrix(d,d,0.6,symm=True) for i,d in enumerate ([random.randint(1,10) for j in range(10)])]
      A = blkdiag(Ai)
      
      #A.sparsity().spy()
      perm =  numpy.random.permutation(range(A.size1()))
      
      AP = A[perm,perm]
      #AP.sparsity().spy()

      n,p,r = AP.sparsity().stronglyConnectedComponents()
      
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
  
    A = DMatrix.sparse(6,4)
    A[0,0] = 1
    A[1,2] = 1
    A[2,2] = 1
    A[5,3] = 1

    ret, rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock = A.sparsity().dulmageMendelsohn()

    # Checked with CSparse
    self.checkarray(DMatrix([ret]),DMatrix([4]))
    self.checkarray(rowperm,DMatrix([2, 3, 4, 1, 0, 5]).T)
    self.checkarray(colperm,DMatrix([ 2,0,3,1]).T)
    self.checkarray(rowblock,DMatrix([ 0, 4,5,6,6]).T)
    self.checkarray(colblock,DMatrix([ 0, 1,2,3,4]).T)
    self.checkarray(coarse_rowblock,DMatrix([ 0, 3,4,6,6]).T)
    self.checkarray(coarse_colblock,DMatrix([ 0, 1,3,3,4]).T)
    
    
    A = DMatrix.sparse(6,4)
    A[0,0] = 1
    A[1,2] = 1
    A[2,2] = 1
    A[5,3] = 1
    A[4,1] = 1
    A[3,0] = 1

    A.sparsity().spy()

    ret, rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock = A.sparsity().dulmageMendelsohn()

    # Checked with CSparse
    self.checkarray(DMatrix([ret]),DMatrix([3]))
    self.checkarray(rowperm,DMatrix([2,3,0,1,4,5]).T)
    self.checkarray(colperm,DMatrix([ 0, 2, 1, 3]).T)
    self.checkarray(rowblock,DMatrix([ 0, 4,5,6]).T)
    self.checkarray(colblock,DMatrix([ 0, 2,3,4]).T)
    self.checkarray(coarse_rowblock,DMatrix([ 0, 2, 4,6,6]).T)
    self.checkarray(coarse_colblock,DMatrix([ 0, 2,4,4,4]).T)
    
    A = DMatrix.sparse(6,4)
    A[0,0] = 1
    A[1,2] = 1
    A[2,2] = 1
    A[5,3] = 1
    A[4,1] = 1
    A[3,0] = 1
    A = A + DMatrix.eye(6)[:,:4]

    A.sparsity().spy()

    ret, rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock = A.sparsity().dulmageMendelsohn()

    # Checked with CSparse
    self.checkarray(DMatrix([ret]),DMatrix([1]))
    self.checkarray(rowperm,DMatrix([4, 5, 0, 1, 2, 3]).T)
    self.checkarray(colperm,DMatrix([ 0, 1, 2, 3]).T)
    self.checkarray(rowblock,DMatrix([ 0, 6]).T)
    self.checkarray(colblock,DMatrix([ 0, 4]).T)
    self.checkarray(coarse_rowblock,DMatrix([ 0, 2, 6,6,6]).T)
    self.checkarray(coarse_colblock,DMatrix([ 0, 4,4,4,4]).T)

    
  def test_jacsparsityHierarchical(self):

    X = SX.sym("X",100)
    P = SX.sym("P",1000)

    optvar = vertcat([X,P])

    p = SX.sym("p")

    g = SXFunction([optvar,p],[X*p])
    g.setOption("verbose",True)
    g.init()

    J = g.jacobian()
    J.setOption("verbose",True)
    J.init()
    
    self.assertTrue(J.getOutput()[:,:X.size()].sparsity()==Sparsity.diag(100))

    X = SX.sym("X",100)
    P = SX.sym("P",1000)

    p = SX.sym("p")

    g = SXFunction([X,p],[vertcat([X*p,P])])
    g.setOption("verbose",True)
    g.init()

    J = g.jacobian()
    J.setOption("verbose",True)
    J.init()
    
    self.assertTrue(J.getOutput()[:X.size(),:].sparsity()==Sparsity.diag(100))
    
  def test_rowcol(self):
    n = 3
    
    s = Sparsity.rowcol([n-1,0],[0,n-1],n,n)
    self.checkarray(IMatrix(s.colind()),IMatrix([0,2,2,4]))
    self.checkarray(IMatrix(s.row()),IMatrix([0,2,0,2]))

  def test_inverse(self):
    numpy.random.seed(0)
    d = self.randDMatrix(20,20,0.6,symm=True)
    sp = d.sparsity()
    
    for sp in [sp,Sparsity.dense(4,4),Sparsity.sparse(4,4),Sparsity.tril(4),Sparsity.tril(4).T]:
    
      d = IMatrix(sp,1)
      
      dt = sparse(1-d)
      dt = IMatrix(dt.sparsity(),1)
      
      trial = IMatrix(sp.patternInverse(),1)
      
      d.printDense()
      dt.printDense()
      trial.printDense()
      
      self.checkarray(trial,dt)
      
  def test_internalapi(self):
    s = Sparsity.dense(2,2)
    with warnings.catch_warnings():
      warnings.simplefilter("error",SyntaxWarning)
      with self.assertRaises(Exception):
        s.reCache()
      
      with internalAPI():
        s.reCache()
    
if __name__ == '__main__':
    unittest.main()

