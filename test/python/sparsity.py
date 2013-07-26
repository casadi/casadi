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
    
    a = CRSSparsity(4,5)
    for i in nza:
      a.getNZ(i[0],i[1])
      
    b = CRSSparsity(4,5)  
    for i in nzb:
      b.getNZ(i[0],i[1])
      
    w = UCharVector()
    c=a.patternUnion(b,w)
    self.assertEquals(w.size(),len(nza.union(nzb)))
    for k in range(w.size()):
      ind = (c.getRow()[k],c.col(k))
      if (ind in nza and ind in nzb):
        self.assertEquals(w[k],1 | 2)
      elif (ind in nza):
        self.assertEquals(w[k],1)
      elif (ind in nzb):
        self.assertEquals(w[k],2)
        
    c = a + b
    self.assertEquals(c.size(),len(nza.union(nzb)))
    for k in range(c.size()):
      ind = (c.getRow()[k],c.col(k))
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
    
    a = CRSSparsity(4,5)
    for i in nza:
      a.getNZ(i[0],i[1])
      
    b = CRSSparsity(4,5)  
    for i in nzb:
      b.getNZ(i[0],i[1])
    
    c=a.patternIntersection(b)
    for k in range(c.size()):
      ind = (c.getRow()[k],c.col(k))
      self.assertTrue(ind in nza and ind in nzb)
        
    c = a * b
    self.assertEquals(c.size(),len(nza.intersection(nzb)))
    for k in range(c.size()):
      ind = (c.getRow()[k],c.col(k))
      self.assertTrue(ind in nza and ind in nzb)
       
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
    
  def tomatrix(self,s):
    d = DMatrix(s,1)
    for k in range(d.size()):
      d[k] = k+1
    return d

  def test_NZ(self):
    self.message("NZ constructor")
    nza = [  (0,0),
             (0,1),
             (2,0),
             (2,3),
             (2,4),
             (3,1)]
    
    a = CRSSparsity(4,5)
    for i in nza:
      a.getNZ(i[0],i[1])
      
    b = sp_triplet(4,5,[i[0] for i in nza],[i[1] for i in nza])
    self.checkarray(self.tomatrix(a),self.tomatrix(b),"rowcol")

  def test_rowcol(self):
    self.message("rowcol constructor")
    
    r = [0,1,3]
    c = [1,4]
    a = CRSSparsity(4,5)
    for i in r:
      for j in c:
        a.getNZ(i,j)
      
    b = sp_rowcol(r,c,4,5)
    self.checkarray(self.tomatrix(a),self.tomatrix(b),"rowcol")
     
  def test_reshape(self):
    self.message("Reshape")
    nza = set([  (0,0),
             (0,1),
             (2,0),
             (2,3),
             (2,4),
             (3,1)])
    
    a = CRSSparsity(4,5)
    for i in nza:
      a.getNZ(i[0],i[1])
      
    A=self.tomatrix(a).toArray()
    B=self.tomatrix(casadi.reshape(a,2,10)).toArray()
    B_=numpy.reshape(A,(2,10))
    
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
    
    a = CRSSparsity(4,5)
    for i in nza:
      a.getNZ(i[0],i[1])
      
    A=self.tomatrix(a).toArray()
    B=self.tomatrix(vec(a)).toArray()
    B_=numpy.reshape(A,(20,1))
    
    self.checkarray(B,B_,"reshape")
    
    
  def test_refcount(self):
      x = DMatrix(sp_tril(4),5)
      s = mul(x,x).sparsity()
      self.assertEqual(s.numel(),16)
      
  def test_splower(self):
    sp = CRSSparsity(3,4,[1,2,1],[0,2,2,3])
    print array(sp)
    print array(lowerSparsity(sp))
    print lowerNZ(sp)
    
    
  def test_diag(self):
    self.message("diag")
    A = CRSSparsity(5,5)
    A.getNZ(1,1)
    A.getNZ(2,4)
    A.getNZ(3,3)
    
    sp, mapping = A.diag()
    B = DMatrix(sp,1)
    
    self.checkarray(array([[0],[1],[0],[1],[0]]),B,"diag(matrix)")
    self.checkarray(array([0,2]),array(list(mapping)),"diag(vector)")
    
    #print B
    
    A = CRSSparsity(5,1)
    A.getNZ(1,0)
    A.getNZ(2,0)
    A.getNZ(4,0)
    
    sp, mapping = A.diag()
    B = DMatrix(sp,1)
    
    self.checkarray(array([[0,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,0,0],[0,0,0,0,1]]),B,"diag(vector)")
    
    self.checkarray(array([0,1,2]),array(list(mapping)),"diag(vector)")
    
    A = CRSSparsity(1,5)
    A.getNZ(0,1)
    A.getNZ(0,2)
    A.getNZ(0,4)
    
    sp, mapping = A.diag()
    B = DMatrix(sp,1)
    
    self.checkarray(array([[0,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,0,0],[0,0,0,0,1]]),B,"diag(vector)")
    
    self.checkarray(array([0,1,2]),array(list(mapping)),"diag(vector)")
    
  def test_vecMX(self):
    self.message("vec MXFunction")
    q = DMatrix([[1,2,3,4,9],[5,6,7,8,8],[9,10,11,12,6],[1,2,3,4,5]])

    X = MX("X",4,5)

    F = MXFunction([X],[X**2])
    F.init()
    F.setInput(q,0)
    F.evaluate()
    F_ = vec(F.getOutput(0))

    G = vec(F)
    G.setInput(vec(q),0)
    G.evaluate()
    G_ = G.getOutput()

    self.checkarray(F_,G_,"vec MX")
    
  def test_vecSX(self):
    self.message("vec SXFunction")
    q = DMatrix([[1,2,3,4,9],[5,6,7,8,8],[9,10,11,12,6],[1,2,3,4,5]])

    X = ssym("X",4,5)

    F = SXFunction([X],[X**2])
    F.init()
    F.setInput(q,0)
    F.evaluate()
    F_ = vec(F.getOutput(0))

    G = vec(F)
    G.setInput(vec(q),0)
    G.evaluate()
    G_ = G.getOutput()

    self.checkarray(F_,G_,"vec SX")
    
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
    
    a = CRSSparsity(5,5)
    for i in nza:
      a.getNZ(i[0],i[1])
      
    b = ssym("b",a)
    
    self.assertRaises(Exception,lambda: b[sp_diag(3)])
    
    d = sp_diag(5)
    c = b[d]

    self.assertTrue(c.sparsity()==d)
    
    f = SXFunction([b],[c])
    f.init()
    f.setInput(range(1,len(nza)+1))
    f.evaluate()
    
    self.checkarray(DMatrix(f.output().data()),DMatrix([1,0,0,7,0]),"sparsity index")
    
    self.assertTrue(f.output().data()[1]==0)
    
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
    
    a = CRSSparsity(5,5)
    for i in nza:
      a.getNZ(i[0],i[1])
      
    b = msym("b",a)
    
    self.assertRaises(Exception,lambda: b[sp_diag(3)])
    
    d = sp_diag(5)
    c = b[d]

    self.assertTrue(c.sparsity()==d)
    
    f = MXFunction([b],[c])
    f.init()
    f.setInput(range(1,len(nza)+1))
    f.evaluate()
    
    self.checkarray(DMatrix(f.output().data()),DMatrix([1,0,0,7,0]),"sparsity index")
    
  def test_getSparsityCCS(self):
    self.message("CCS format")
    nza = set([  (0,0),
             (0,1),
             (2,0),
             (2,3),
             (3,3),
             (2,4),
             (3,1)])
    
    a = CRSSparsity(4,5)
    for i in nza:
      a.getNZ(i[0],i[1])
      
    A1, B1= a.getSparsityCCS()
    
    A2, B2 = (a.T).getSparsityCRS()
    
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
      
      ST = Ar.T.sparsity()
      
      blocks = []
      acc = -1
      mc = 0
      for i in range(0,Ar.size1()):
        mc = max(ST.col()[ST.rowind()[i+1]-1],mc)
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
       
      ST = Ar.T.sparsity()
      
      blocks = []
      acc = -1
      mc = 0
      for i in range(0,Ar.size1()):
        mc = max(ST.col()[ST.rowind()[i+1]-1],mc)
        if mc==i:
          blocks.append(i-acc)
          acc = i
      
      truth = [i.size1() for i in Ai]
      tryme = blocks

      self.assertTrue(n>=len(truth))
      self.assertTrue(n>=len(tryme))
    
  def test_dm(self):
  
    A = DMatrix(6,4)
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
    
    
    A = DMatrix(6,4)
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
    
    A = DMatrix(6,4)
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

    X = ssym("X",100)
    P = ssym("P",1000)

    optvar = vertcat([X,P])

    p = ssym("p")

    g = SXFunction([optvar,p],[X*p])
    g.setOption("verbose",True)
    g.init()

    J = g.jacobian()
    J.setOption("verbose",True)
    J.init()
    
    self.assertTrue(J.output()[:,:X.size()].sparsity()==sp_diag(100))

    X = ssym("X",100)
    P = ssym("P",1000)

    p = ssym("p")

    g = SXFunction([X,p],[vertcat([X*p,P])])
    g.setOption("verbose",True)
    g.init()

    J = g.jacobian()
    J.setOption("verbose",True)
    J.init()
    
    self.assertTrue(J.output()[:X.size(),:].sparsity()==sp_diag(100))


if __name__ == '__main__':
    unittest.main()

