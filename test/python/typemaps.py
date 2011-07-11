from casadi import *
import casadi as c
from numpy import *
import numpy as n
import unittest
from types import *
from helpers import *
from scipy.sparse import *
from scipy import *
import warnings

class typemaptests(casadiTestCase):

  def setUp(self):
    pass

  def test_0(self):
    self.message("Typemap array -> DMatrix")
    arrays = [array([[1,2,3],[4,5,6]]),array([[1,2],[3,4],[5,6]],dtype=double),array([[3.2,4.6,9.9]])]
    for i in range(len(arrays)):
      m = arrays[i]
      zt=trans(trans(m))
      self.assertTrue(isinstance(zt,DMatrix),"DMatrix expected")
      self.checkarray(m,zt,"DMatrix(numpy.ndarray)")
      self.checkarray(m,zt.toArray(),"DMatrix(numpy.ndarray).toArray()")
      self.checkarray(m,zt.toCsr_matrix(),"DMatrix(numpy.ndarray).toCsr_matrix()")
      
  def test_1(self):
    self.message("DMatrix -> DMatrix")
    arrays = [DMatrix(3,4,[1,2,1],[0,2,2,3],[3,2.3,8])]
    for i in range(len(arrays)):
      m = arrays[i]
      zt=trans(trans(m))
      self.assertTrue(isinstance(zt,DMatrix),"DMatrix expected")
      self.checkarray(m,zt,"DMatrix(DMatrix)")
      self.checkarray(m,zt.toArray(),"DMatrix(DMatrix).toArray()")
      self.checkarray(m,zt.toCsr_matrix(),"DMatrix(DMatrix).toCsr_matrix()")
   
  def test_2(self):
    self.message("crs_matrix -> DMatrix")
    arrays = [csr_matrix( ([3,2.3,8],([0,2,0],[1,1,2])), shape = (3,4), dtype=double ),
              csr_matrix( ([3,2.3,8],([0,2,0],[1,1,2])), shape = (3,4), dtype=int )
              ]
    for i in range(len(arrays)):
      m = arrays[i]
      zt=trans(trans(m))
      self.assertTrue(isinstance(zt,DMatrix),"DMatrix expected")
      self.checkarray(m,zt,"DMatrix(crs_matrix)")
      self.checkarray(m,zt.toArray(),"DMatrix(crs_matrix).toArray()")
      self.checkarray(m,zt.toCsr_matrix(),"DMatrix(crs_matrix).toCsr_matrix()")
      
      
  def test_setget(self):
    self.message("DMatrix set/get")
    data = n.array([3,2.3,8])
    dm=DMatrix(3,4,[1,2,1],[0,2,2,3],[3,2.3,8])
    
    c=dm.toCsr_matrix()
    z=n.zeros((3,4))
    dm.get(z)
    self.checkarray(z,dm,"get(2Dndarray)")
    z=n.matrix(n.zeros((3,4)))
    dm.get(z)
    self.checkarray(z,dm,"get(2Dmatrix)")
    z=n.zeros((12,5))
    self.assertRaises(TypeError,lambda : dm.get(z),"get(wrong size ndarray)")
    z=ones((3,4))
    dm.set(z)
    self.checkarray(dm.toArray() > 0,dm,"set(2Dndarray)")
    z=n.matrix(ones((3,4)))
    dm.set(z)
    self.checkarray(dm.toArray() > 0,dm,"set(2Dmatrix)")
    z=n.zeros((12,5))
    self.assertRaises(TypeError,lambda : dm.set(z))
    
    dm.set(c)
    self.checkarray(c,dm,"set(csr_matrix)")
    
    z=n.zeros(3)
    dm.get(z)
    self.checkarray(n.matrix(z),n.matrix(data),"get(1Dndarray)")
    dm.set(z)
    self.checkarray(c,dm,"set(1Dndarray)")

    dm = dm * 2
    dm.get(c)
    dm.shape = (dm.size1(),dm.size2())

    self.checkarray(c,dm,"get(csr_matrix)")
    
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      c[0,0]=1
    self.assertRaises(TypeError,lambda :  dm.set(c))
    self.assertRaises(TypeError,lambda :  dm.get(c))

  def test_conversion(self):
    self.message("DMatrix conversions")
    w = DMatrix(3,4,[1,2,1],[0,2,2,3],[3,2.3,8])
    d = array([[1,2,3],[4,5,6]])
    
    list(w.data())
    tuple(w.data())
    w.toArray()
    array(w)
    w.toMatrix()
    matrix(w)
    w.toCsr_matrix()

    self.checkarray(DMatrix(d),d,"DMatrix(numpy.ndarray)")
    #self.checkarray(DMatrix(array([1,2,3,4,5,6])),d.ravel(),"DMatrix(numpy.ndarray)")
    #print DMatrix(array([1,2,3,4,5,6]))
    #print DMatrix(array([1,2,3,6]),2,2).toArray()

    #print DMatrix(array([1,2,3,6]),2,2).toArray()
    
  def test_autoconversion(self):
    self.message("Auto conversion DMatrix")
    x=array([2.3])
    s = DMatrix([[1,2],[3,4]])
    n = array(s)
    
    self.checkarray(x[0]*s,s*x[0],"")
    self.checkarray(x[0]*s,n*x[0],"")
    
    self.checkarray(x[0]/s,1/(s/x[0]),"")
    self.checkarray(x[0]/s,x[0]/n,"")
    
    self.checkarray(x[0]-s,-(s-x[0]),"")
    self.checkarray(x[0]-s,x[0]-n,"")
    
    self.checkarray(x[0]+s,s+x[0],"")
    self.checkarray(x[0]+s,x[0]+n,"")
    
    w=array([2.3])[0]
    w+=s
    self.checkarray(w,2.3+n,"")

    w=array([2.3])[0]
    w-=s
    self.checkarray(w,2.3-n,"")
    
    w=array([2.3])[0]
    w*=s
    self.checkarray(w,2.3*n,"")
    
    w=array([2.3])[0]
    w/=s
    self.checkarray(w,2.3/n,"")
    
    x=[2.3]

    self.checkarray(x[0]*s,s*x[0],"")
    self.checkarray(x[0]*s,n*x[0],"")
    
    self.checkarray(x[0]/s,1/(s/x[0]),"")
    self.checkarray(x[0]/s,x[0]/n,"")
    
    self.checkarray(x[0]-s,-(s-x[0]),"")
    self.checkarray(x[0]-s,x[0]-n,"")
    
    self.checkarray(x[0]+s,s+x[0],"")
    self.checkarray(x[0]+s,x[0]+n,"")
    
    
    w=2.3
    w+=s
    self.checkarray(w,2.3+n,"")
    
    w=2.3
    w-=s
    self.checkarray(w,2.3-n,"")
    
    w=2.3
    w*=s
    self.checkarray(w,2.3*n,"")
    
    w=2.3
    w/=s
    self.checkarray(w,2.3/n,"")

  def test_autoconversionMX(self):
    self.message("Auto conversion MX")
    s = DMatrix([[1,2],[3,4]])
    x = SX(3)
    y = MX(3)
    
    for (s,x,y) in [
                  (DMatrix([[1,2],[3,4]]),SX("x"),MX("x",1,1)),
                  (3,symbolic("x",2,2),MX("x",2,2)),
                  (DMatrix(3),symbolic("x",2,2),MX("y",2,2)),
                  (array([[1,2],[3,4]]),SX("x"),MX("x",1,1))
                  ]:
      for z in [x,y]:
        print "z = %s, s = %s" % (str(z),str(s))
        print "  z = %s, s = %s" % (type(z),type(s))
        -z
        -s
        z+s
        s+z
        s*z
        z*s
        s-z
        z-s
        z/s
        s/z
        s**z
        z**s
      
  def test_set(self):
    self.message("DMatrix set on dense matrices")
    
    # should be integers
    goallist = [1,2,3]
    goal = array([goallist]).T

    test={
      "list" : goallist,
      "tuple" : tuple(goallist),
      "array1ddouble" : array(goallist,dtype=double),
      "array2ddouble" : array([goallist],dtype=double).T,
      "array1dint" : array(goallist),
      "array2dint" : array([goallist]).T
    }
    w=DMatrix(goal)
    self.checkarray(w,goal,"Constructor")
    
    for name, value in test.items():
      w.set(value)
      self.checkarray(w,goal,"name")
      
      
  def testGenericType(self):
    self.message("Generic type")
    x=SX("x")
    f=SXFunction([x],[2*x])
    f.setOption("name","foo")
    self.assertEquals(f.getOption("name"),"foo")
    f.setOption("verbose",True)
    #self.assertTrue(isinstance(f.getOption("verbose"),bool))
    self.assertTrue(f.getOption("verbose"))
    f.setOption("verbose",False)
    self.assertTrue(not(f.getOption("verbose")))
    f.setOption("number_of_adj_dir",3)
    self.assertTrue(isinstance(f.getOption("number_of_adj_dir"),int))
    self.assertEquals(f.getOption("number_of_adj_dir"),3)
    d=f.dictionary()
    self.assertTrue(isinstance(d,dict))
    d["verbose"]=True
    d["number_of_adj_dir"]=7
    f.setOption(d)
    self.assertTrue(f.getOption("verbose"))
    self.assertEquals(f.getOption("number_of_adj_dir"),7)

  def testGenericType2(self):
    self.message("Generic type 2")
    for i in [0,1,7,-7]:
	    a=GenericType(i)
	    self.assertTrue(a.isInt())
	    self.assertFalse(a.isBool())
	    self.assertFalse(a.isDouble())
	    self.assertFalse(a.isString())
	    self.assertEqual(a.toInt(),i)
    for i in [True,False]:
	    a=GenericType(i)
	    #self.assertFalse(a.isInt())
	    #self.assertTrue(a.isBool())
	    #self.assertFalse(a.isDouble())
	    #self.assertEqual(a.toBool(),i)

    for i in [0.01,-5.7]:
	    a=GenericType(i)
	    self.assertFalse(a.isInt())
	    self.assertFalse(a.isBool())
	    self.assertTrue(a.isDouble())
	    self.assertFalse(a.isString())
	    self.assertEqual(a.toDouble(),i)

    for i in ["","foo"]:
	    a=GenericType(i)
	    self.assertFalse(a.isInt())
	    self.assertFalse(a.isBool())
	    self.assertFalse(a.isDouble())
	    self.assertTrue(a.isString())
	    self.assertEqual(a.toString(),i)

    for i in [(0,1,5)]:
	    a=GenericType(i)
	    #print a
	    #self.assertTrue(a.isIntVector())
	    #self.assertTrue(a.isDoubleVector())
	    #self.assertEqual(a.toString(),i)

  def testGenericType3(self):
    self.message("Generic type 3")
    
    is_differential_ivec = IVector(2)
    is_differential_ivec[0] = 4

    is_differential_gentype = GenericType(is_differential_ivec)
    
    self.assertTrue(is_differential_gentype.isIntVector())

	    
  def test_operators(self):
    self.message("Test operators on mixed numpy.array/Matrix")
    self.message(":SXMatrix")
    x=SX("x")
    y=SX("y")

    C=SXMatrix([x,y])
    N=matrix([x,y]).T
    
    self.assertTrue(isinstance(N+C,SXMatrix))
    self.assertTrue(isinstance(C+N,SXMatrix))
    
    f=SXFunction([[x,y]],[C+N])
    f.init()
    f.input().set([7,13])
    f.evaluate()
    self.checkarray(f.output(),matrix([14,26]).T,"addition")
    
    f=SXFunction([[x,y]],[N+C])
    f.init()
    f.input().set([7,13])
    f.evaluate()
    self.checkarray(f.output(),matrix([14,26]).T,"addition")
    
    self.message(":DMatrix")
    D=DMatrix([7,13])
    N=matrix([7,13]).T
    
    self.assertTrue(isinstance(N+D,DMatrix))
    self.assertTrue(isinstance(D+N,DMatrix))
    
    self.checkarray(N+D,matrix([14,26]).T,"addition")
    self.checkarray(D+N,matrix([14,26]).T,"addition")
    
    self.assertTrue(isinstance(C+D,SXMatrix))
    self.assertTrue(isinstance(D+C,SXMatrix))
  
      
    f=SXFunction([[x,y]],[C+D])
    f.init()
    f.input().set([1,4])
    f.evaluate()
    self.checkarray(f.output(),matrix([8,17]).T,"addition")
    
    f=SXFunction([[x,y]],[D+C])
    f.init()
    f.input().set([1,4])
    f.evaluate()
    self.checkarray(f.output(),matrix([8,17]).T,"addition")

  def test_DMatrixSXMatrixcast(self):
    self.message("Casting DMatrix to SXMatrix")
    
    W = SXMatrix(DMatrix([[1,2,3],[4,5,6]]))

    self.assertEqual(W.size1(),2)
    self.assertEqual(W.size2(),3)

  def test_DMatrixMXcast(self):
    self.message("Casting DMatrix to MX")
    W = MX(DMatrix([[1,2,3],[4,5,6]]))
    
    self.assertEqual(W.size1(),2)
    self.assertEqual(W.size2(),3)
    
  def test_DMatrixSXMatrix(self):
    self.message("Casting DMatrix to SXMatrix")
    
    w = DMatrix([[1,2,3],[4,5,6]])
    x = SX("x")
    
    f = SXFunction([x],[w])
    
    W = f.outputSX()
    self.assertEqual(W.size1(),2)
    self.assertEqual(W.size2(),3)

  def test_DMatrixMX(self):
    self.message("Casting DMatrix to MX")
    w = DMatrix([[1,2,3],[4,5,6]])
    x = MX("x")
    
    f = MXFunction([x],[w])
    
    W = f.outputMX()

    self.assertEqual(W.size1(),2)
    self.assertEqual(W.size2(),3)

  def test_sharedarray(self):
    w = DMatrix([[1,2],[3,4]])
    W = w.toArray(shared=True)
    self.checkarray(w,W,"shared")
    
    w[0,1] = 8
    self.checkarray(w,W,"shared")

    W[:,0] = 47
    self.checkarray(w,W,"shared")
    
  def test_setgetslicetransp(self):
    self.message("set/get on DMatrix using tranpose")
    
    w = DMatrix([[0,0],[0,0]])

    A = matrix([[1.0,2],[3,4]])
    B = matrix([[4.0,5],[6,7]])
    
    w.set(A)
    w.get(B.T)
    
    self.checkarray(B.T,A,"get")
    
  def test_setgetslice(self):
    self.message("set/get on DMatrix using slices")
    
    w = DMatrix([[0,0]])

    A = matrix([[1.0,2],[3,4]])
    B = matrix([[4.0,5],[6,7]])
    
    w.set(A[0,:])
    self.checkarray(w,A[0,:],"set")
    w.get(B[0,:])
    self.checkarray(B[0,:],A[0,:],"get")
    
    w = DMatrix([[0],[0]])


    w.set(A[:,0])
    self.checkarray(w,A[:,0],"set")
    w.get(B[:,0])
    self.checkarray(B[:,0],A[:,0],"get")
    
    w = DMatrix([[1,2],[3,4]])
    A = zeros((8,7))
    B = zeros((8,7))
    w.get(A[2:7:3,:7:4])
    B[2:7:3,:7:4] = w
    self.checkarray(A,B,"get")
    
  def test_vertcatprecedence(self):
    self.message("Argument precedence DMatrix")
    a = DMatrix([1,2])
    self.assertTrue(isinstance(vertcat([a,a]),DMatrix))
    
    a = DMatrix([1,2])
    self.assertTrue(isinstance(vertcat([a,[1,2,3]]),DMatrix))
    
    
    a = MX([1,2])
    print vertcat([a,[1,2,3]])
    self.assertTrue(isinstance(vertcat([a,[1,2,3]]),MX))
    
if __name__ == '__main__':
    unittest.main()
