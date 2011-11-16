from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *
import numpy 

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
      
    #w = IVector()
    #c=a.patternUnion(b,w)
    #self.assertEquals(w.size(),len(nza.union(nzb)))
    #for k in range(w.size()):
      #ind = (c.getRow()[k],c.col(k))
      #if (ind in nza and ind in nzb):
        #self.assertEquals(w[k],1 | 2)
      #elif (ind in nza):
        #self.assertEquals(w[k],1)
      #elif (ind in nzb):
        #self.assertEquals(w[k],2)
        
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
      
    b = sp_NZ([i[0] for i in nza],[i[1] for i in nza],4,5,True)
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
      return #Ticket 147
      x = DMatrix(sp_tril(4),5)
      s = c.mul(x,x).sparsity()
      self.assertEqual(s.numel(),10)
      
  def test_splower(self):
    sp = CRSSparsity(3,4,[1,2,1],[0,2,2,3])
    print array(sp)
    print array(lowerSparsity(sp))
    print lowerNZ(sp)
    
    
  def test_diag(self):
    self.message("diag")
    mapping = IVector()
    A = CRSSparsity(5,5)
    A.getNZ(1,1)
    A.getNZ(2,4)
    A.getNZ(3,3)
    
    B = DMatrix(A.diag(mapping),1)
    
    self.checkarray(array([[0],[1],[0],[1],[0]]),B,"diag(matrix)")
    self.checkarray(array([0,2]),array(list(mapping)),"diag(vector)")
    
    #print B
    
    A = CRSSparsity(5,1)
    A.getNZ(1,0)
    A.getNZ(2,0)
    A.getNZ(4,0)
    
    B = DMatrix(A.diag(mapping),1)
    
    self.checkarray(array([[0,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,0,0],[0,0,0,0,1]]),B,"diag(vector)")
    
    self.checkarray(array([0,1,2]),array(list(mapping)),"diag(vector)")
    
    A = CRSSparsity(1,5)
    A.getNZ(0,1)
    A.getNZ(0,2)
    A.getNZ(0,4)
    
    B = DMatrix(A.diag(mapping),1)
    
    self.checkarray(array([[0,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,0,0],[0,0,0,0,1]]),B,"diag(vector)")
    
    self.checkarray(array([0,1,2]),array(list(mapping)),"diag(vector)")
    
      
if __name__ == '__main__':
    unittest.main()

