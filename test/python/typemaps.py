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
import numpy as n
import unittest
from types import *
from helpers import *

scipy_available = True
try:
	from scipy.sparse import *
	from scipy import *
except:
	scipy_available = False
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
      if scipy_available:
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
      if scipy_available:
        self.checkarray(m,zt.toCsr_matrix(),"DMatrix(DMatrix).toCsr_matrix()")
   
  def test_2(self):
    self.message("crs_matrix -> DMatrix")
    if not(scipy_available):
      return
    arrays = [csr_matrix( ([3,2.3,8],([0,2,0],[1,1,2])), shape = (3,4), dtype=double ),
              csr_matrix( ([3,2.3,8],([0,2,0],[1,1,2])), shape = (3,4), dtype=int )
              ]
    for i in range(len(arrays)):
      m = arrays[i]
      zt=trans(trans(m))
      self.assertTrue(isinstance(zt,DMatrix),"DMatrix expected")
      self.checkarray(m,zt,"DMatrix(crs_matrix)")
      self.checkarray(m,zt.toArray(),"DMatrix(crs_matrix).toArray()")
      if scipy_available:
        self.checkarray(m,zt.toCsr_matrix(),"DMatrix(crs_matrix).toCsr_matrix()")
      
      
  def test_setget(self):
    self.message("DMatrix set/get")
    data = n.array([3,2.3,8])
    dm=DMatrix(3,4,[1,2,1],[0,2,2,3],[3,2.3,8])
    
    if scipy_available:
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
    
    if scipy_available:
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
    if scipy_available:
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
    
    def doit(z,s,fun):
      function = None
      
      if type(z) in [type(SX()),type(SXMatrix())]:
        ztype = [type(SX()),type(SXMatrix())]
        function = SXFunction
      
      if type(z) in [type(MX())]:
        ztype = [type(MX())]
        function = MXFunction
        
      r = fun(z,s)
            
      if type(z) is type(SX()) and type(s) is type(SX()):
        self.assertTrue(type(r) is type(SX()))
        

      self.assertTrue(type(r) in ztype,"Expected %s but got %s" % (str(ztype),str(type(r))))
      
      hasNum = True
      if type(s) in [type(SX()),type(MX()),type(SXMatrix())]:
        hasNum = False
      
      if hasNum:
        dummy = [1.3,2.7,9.4,1.0]

        f=function([z],[r])
        f.init()
        f.input().set(dummy[0:f.input().size()])
        f.evaluate()
        
        f_=function([z],[z])
        f_.init()
        f_.input().set(dummy[0:f.input().size()])
        f_.evaluate()
        

        self.checkarray(fun(f_.output(),DMatrix(s)),f.output(),"operation")
      else:
        dummy = [1.3,2.7,9.4,1.0]
        dummy2 = [0.3,2.4,1.4,1.7]
        
        f=function([z,s],[r])
        f.init()
        f.input(0).set(dummy[0:f.input(0).size()])
        f.input(1).set(dummy2[0:f.input(1).size()])
        f.evaluate()
        
        f_=function([z,s],[z,s])
        f_.init()
        f_.input(0).set(dummy[0:f.input(0).size()])
        f_.input(1).set(dummy2[0:f.input(1).size()])
        f_.evaluate()

        self.checkarray(fun(f_.output(0),f_.output(1)),f.output(),"operation")
    
    
    def tests(z,s):
      doit(z,s,lambda z,s: -z)
      doit(z,s,lambda z,s: z+s)
      doit(z,s,lambda z,s: s+z)
      doit(z,s,lambda z,s: s*z)
      doit(z,s,lambda z,s: z*s)
      doit(z,s,lambda z,s: z-s)
      doit(z,s,lambda z,s: s-z)
      doit(z,s,lambda z,s: z/s)
      doit(z,s,lambda z,s: s/z)
      doit(z,s,lambda z,s: z**s)
      doit(z,s,lambda z,s: s**z)
      doit(z,s,lambda z,s: fmin(s,z))
      doit(z,s,lambda z,s: fmax(s,z))
      doit(z,s,lambda z,s: min(s,z))
      doit(z,s,lambda z,s: max(s,z))
      doit(z,s,lambda z,s: constpow(s,z))
      doit(z,s,lambda z,s: constpow(z,s))
      doit(z,s,lambda z,s: arctan2(s,z))
      doit(z,s,lambda z,s: arctan2(z,s))

    nums = [array([[1,2],[3,4]]),DMatrix([[1,2],[3,4]]), DMatrix(4), array(4),4.0,4]
        
    ## numeric & SXMatrix
    for s in nums:
      for z in [SX("x"), ssym("x"), ssym("x",2,2)]:
        print "z = %s, s = %s" % (str(z),str(s))
        print "  z = %s, s = %s" % (type(z),type(s))
        tests(z,s)
       
    # numeric & MX
    for s in nums:
      for z in [MX("x",2,2)]:
        print "z = %s, s = %s" % (str(z),str(s))
        print "  z = %s, s = %s" % (type(z),type(s))
        tests(z,s)
        
    # SX & SX
    for s in [SX("x"), ssym("x"), ssym("x",2,2)]:
      for z in [SX("x"),ssym("x"), ssym("x",2,2)]:
        print "z = %s, s = %s" % (str(z),str(s))
        print "  z = %s, s = %s" % (type(z),type(s))
        tests(z,s)
         
    ## MX & MX
    for s in [MX("x"),MX("x",2,2)]:
      for z in [MX("x"),MX("x",2,2)]:
        print "z = %s, s = %s" % (str(z),str(s))
        print "  z = %s, s = %s" % (type(z),type(s))
        tests(z,s)
        
    for (s,x,y) in [
                  (matrix([[1,2],[3,4]]),ssym("x",2,2),MX("x",2,2))    
                  ]:
      for z,ztype in zip([x,y],[[type(SXMatrix()),type(SX())],[type(MX())]]):
        print "z = %s, s = %s" % (str(z),str(s))
        print "  z = %s, s = %s" % (type(z),type(s))
        doit(z,s,lambda z,s: -z)
        -s
        doit(z,s,lambda z,s: z+s)
        doit(z,s,lambda s,z: s+z)
        doit(z,s,lambda z,s: s*z)
        doit(z,s,lambda s,z: z*s)
        doit(z,s,lambda z,s: z-s)
        doit(z,s,lambda s,z: s-z)
        doit(z,s,lambda z,s: z/s)
        doit(z,s,lambda s,z: s/z)
        
  def test_conversion_operators(self):
    self.message("COnversion operations")
    

    def doit(z,s,fun):
      function = None
      
      if type(z) in [type(SX()),type(SXMatrix())]:
        ztype = [type(SX()),type(SXMatrix())]
        function = SXFunction
      
      if type(z) in [type(MX())]:
        ztype = [type(MX())]
        function = MXFunction
        
      r = fun(z,s)
            
      if type(z) is type(SX()) and type(s) is type(SX()):
        self.assertTrue(type(r) is type(SX()))
        

      self.assertTrue(type(r) in ztype,"Expected %s but got %s" % (str(ztype),str(type(r))))
      
    def tests(z,s):
      doit(z,s,lambda z,s: s>=z)
      doit(z,s,lambda z,s: s>z)
      doit(z,s,lambda z,s: s<=z)
      doit(z,s,lambda z,s: s<z)
      doit(z,s,lambda z,s: s==z)
      doit(z,s,lambda z,s: s!=z)
      
    nums = [array([[1,2],[3,4]]),DMatrix([[1,2],[3,4]]), DMatrix(4), array(4),4.0,4]
        
    ## numeric & SXMatrix
    for s in nums:
      for z in [SX("x"), ssym("x"), ssym("x",2,2)]:
        print "z = %s, s = %s" % (str(z),str(s))
        print "  z = %s, s = %s" % (type(z),type(s))
        tests(z,s)
       
    # numeric & MX
    for s in nums:
      for z in [MX("x",2,2)]:
        print "z = %s, s = %s" % (str(z),str(s))
        print "  z = %s, s = %s" % (type(z),type(s))
        tests(z,s)
        
    # SX & SX
    for s in [SX("x"), ssym("x"), ssym("x",2,2)]:
      for z in [SX("x"),ssym("x"), ssym("x",2,2)]:
        print "z = %s, s = %s" % (str(z),str(s))
        print "  z = %s, s = %s" % (type(z),type(s))
        tests(z,s)
         
    # MX & MX
    for s in [MX("x"),MX("x",2,2)]:
      for z in [MX("x"),MX("x",2,2)]:
        print "z = %s, s = %s" % (str(z),str(s))
        print "  z = %s, s = %s" % (type(z),type(s))
        tests(z,s)
        
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
      "array2dint" : array([goallist]).T,
      "mixed" : [1,DMatrix(2),array(3)]
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
	    self.assertTrue(a.isIntVector())
	    self.assertFalse(a.isDoubleVector())

    for i in [(0.3,1,5)]:
	    a=GenericType(i)
	    self.assertFalse(a.isIntVector())
	    self.assertTrue(a.isDoubleVector())
	    
    a = GenericType(["foo","bar"])
    self.assertTrue(a.isStringVector())
    x = SX("x")
    f = SXFunction([x],[x])
    #f.setOption("monitor",["foo","bar"])
    #self.assertEqual(f.getOption("monitor")[0],"foo")
    #self.assertEqual(f.getOption("monitor")[1],"bar")
    #f.setOption("monitor",["foo"])
    #self.assertEqual(f.getOption("monitor")[0],"foo")
    #f.setOption("monitor",[])
    
    t=SX("t")

    x=SX("x") 
    dx=SX("dx")

    f=SXFunction(daeIn(t=t, x=vertcat([x,dx])),[vertcat([dx,-x])])
    f.init()

    integrator = CVodesIntegrator(f)
    integrator.setOption("fsens_scaling_factors",[5.0,7])
    integrator.setOption("fsens_scaling_factors",[])
    
  def testGenericType3(self):
    self.message("Generic type 3")
    
    is_differential_ivec = IVector(2)
    is_differential_ivec[0] = 4

    is_differential_gentype = GenericType(is_differential_ivec)
    
    self.assertTrue(is_differential_gentype.isIntVector())

  @requires("IpoptSolver")
  def testGenericTypeBoolean(self):
    x=SX("x")

    nlp = SXFunction(nlIn(x=x),nlOut(f=x**2))
    nlp.init()

    nlp_solver = IpoptSolver(nlp)
    
    self.assertRaises(RuntimeError,lambda : nlp_solver.setOption('acceptable_tol',SX("x")))
    nlp_solver.setOption('acceptable_tol',DMatrix(1))
	    
  def test_operators(self):
    self.message("Test operators on mixed numpy.array/Matrix")
    self.message(":SXMatrix")
    x=SX("x")
    y=SX("y")

    C=SXMatrix([x,y])
    N=matrix([x,y]).T
    
    self.assertTrue(isinstance(N+C,SXMatrix))
    self.assertTrue(isinstance(C+N,SXMatrix))
    
    f=SXFunction([vertcat([x,y])],[C+N])
    f.init()
    f.input().set([7,13])
    f.evaluate()
    self.checkarray(f.output(),matrix([14,26]).T,"addition")
    
    f=SXFunction([vertcat([x,y])],[N+C])
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
  
      
    f=SXFunction([vertcat([x,y])],[C+D])
    f.init()
    f.input().set([1,4])
    f.evaluate()
    self.checkarray(f.output(),matrix([8,17]).T,"addition")
    
    f=SXFunction([vertcat([x,y])],[D+C])
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
    
    W = f.outputExpr(0)
    self.assertEqual(W.size1(),2)
    self.assertEqual(W.size2(),3)

  def test_DMatrixMX(self):
    self.message("Casting DMatrix to MX")
    w = DMatrix([[1,2,3],[4,5,6]])
    x = MX("x")
    
    f = MXFunction([x],[w])
    
    W = f.outputExpr(0)

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
    
  def test_issue190(self):
    self.message("regression test issue #190")
    x=SX("x")
    x * numpy.array(1)
    x * numpy.array(1.2)

    ssym("x") * numpy.array(1.0) 
    MX("x") * numpy.array(1.0)
    
  def test_array_cat(self):
    horzcat((ssym("x",4,3),ones((4,3))))
    
    
  def test_issue(self):
    self.message("std::vector<double> typemap.")
    a = array([0,2,2,3])
    b = array([0.738,0.39,0.99])
    DMatrix(3,4,[1,2,1],[0,2,2,3],[0.738,0.39,0.99])
    DMatrix(3,4,[1,2,1],(0,2,2,3),[0.738,0.39,0.99])
    DMatrix(3,4,[1,2,1],list(a),[0.738,0.39,0.99])
    DMatrix(3,4,[1,2,1],a,[0.738,0.39,0.99])
    DMatrix(3,4,[1,2,1],[0,2,2,3],(0.738,0.39,0.99))
    DMatrix(3,4,[1,2,1],[0,2,2,3],list(b))
    DMatrix(3,4,[1,2,1],[0,2,2,3],b)
    
  def test_imatrix(self):
    self.message("IMatrix")
    
    A = IMatrix(2,2,1)
    B = A + 1
    self.assertEqual(type(B),type(A))
    self.checkarray(array(B),DMatrix(2,2,2),"Imatrix")
    
  def test_issue314(self):
    self.message("regression test for #314: SXMatrix sparsity constructor")
    SXMatrix(sp_diag(3),[1,2,3])
  def test_setAll_365(self):
    self.message("ticket #365: DMAtrix.setAll does not work for 1x1 Matrices as input")
    m = DMatrix.ones(5,5)
    m.setAll(DMatrix(4))
    m.setAll(IMatrix(4))
    
  def test_issue_(self):
    self.message("Ticket #533")


    x=ssym("x")
    z=ssym("z")
    p=ssym("p")
    f = SXFunction(daeIn(x=x,p=p,z=z),daeOut(ode=z/p,alg=z-x))
    f.init()

    integr = IdasIntegrator(f)
    integr.setOption("init_xdot",GenericType())
    integr.setOption("calc_icB",True)
    integr.setOption("t0",0.2)
    integr.setOption("tf",2.3)
    integr.init()

    integr.input("x0").set(7.1)
    integr.input("p").set(2)

    integr.evaluate(1,1)

    print integr.dictionary()


    integr = IdasIntegrator(f)
    integr.setOption("init_xdot",None)
    integr.setOption("calc_icB",True)
    integr.setOption("t0",0.2)
    integr.setOption("tf",2.3)
    integr.init()

    integr.input("x0").set(7.1)
    integr.input("p").set(2)

    integr.evaluate(1,1)

    print integr.dictionary()

    integr = IdasIntegrator(f)
    integr.setOption("init_xdot",[7.1])
    integr.setOption("calc_icB",True)
    integr.setOption("augmented_options", {"init_xdot":GenericType()})

    integr.setOption("t0",0.2)
    integr.setOption("tf",2.3)
    integr.init()

    integr.input("x0").set(7.1)
    integr.input("p").set(2)

    integr.evaluate(1,1)

    print integr.dictionary()

    integr = IdasIntegrator(f)
    integr.setOption("init_xdot",[7.1])
    integr.setOption("calc_icB",True)
    integr.setOption("augmented_options", {"init_xdot":None})

    integr.setOption("t0",0.2)
    integr.setOption("tf",2.3)
    integr.init()

    integr.input("x0").set(7.1)
    integr.input("p").set(2)

    integr.evaluate(1,1)

    print integr.dictionary()
    
  def test_issue570(self):
    self.message("Issue #570: long int")
    longint = 10**50
    print type(longint)
    print casadi.SX('x') + longint
    print longint + casadi.SX('x')
    print casadi.ssym('x') + longint
    print longint + casadi.ssym('x')
    
  def test_casting_DMatrix(self):
    self.message("casting DMatrix")
    
    x = ssym("x")
    f = SXFunction([x],[x])
    f.init()
    class Foo:
      def __DMatrix__(self):
        return DMatrix([4])
        
    f.setInput(Foo())
    self.assertEqual(f.input(),4)

    class Foo:
      def __DMatrix__(self):
        return SXMatrix([4])
        
    self.assertRaises(TypeError,lambda :f.setInput(Foo()))
    
    class Foo:
      def __DMatrix__(self):
        raise Exception("15")
        
    self.assertRaises(TypeError,lambda :f.setInput(Foo()))

    class Foo:
      pass
        
    self.assertRaises(NotImplementedError,lambda :f.setInput(Foo()))

  def test_casting_IMatrix(self):
    self.message("casting IMatrix")

    class Foo:
      def __IMatrix__(self):
        return IMatrix([[4,6],[2,4]])
        
    self.assertEqual(det(Foo()),4)

    class Foo:
      def __IMatrix__(self):
        return SXMatrix([[4,6],[2,4]])
        
    self.assertRaises(TypeError,lambda :det(Foo()))
    
    class Foo:
      def __IMatrix__(self):
        raise Exception("15")
        
    self.assertRaises(TypeError,lambda :det(Foo()))

    class Foo:
      pass
        
    self.assertRaises(NotImplementedError,lambda : det(Foo()))

  def test_casting_SXMatrix(self):
    self.message("casting SXMatrix")
    
    
    x = ssym("x")
    
    class Foo:
      def __SXMatrix__(self):
        return x
        
    SXFunction([x],[Foo()])
    
    class Foo:
      def __SXMatrix__(self):
        return MX("x")
        
    self.assertRaises(TypeError,lambda : SXFunction([x],[Foo()]))
    
    class Foo:
      def __SXMatrix__(self):
        raise Exception("15")
        
    self.assertRaises(TypeError,lambda : SXFunction([x],[Foo()]))

    class Foo:
      pass
        
    self.assertRaises(NotImplementedError,lambda :SXFunction([x],[Foo()]))


  def test_casting_MX(self):
    self.message("casting MX")
    
    
    x = msym("x")
    
    class Foo:
      def __MX__(self):
        return x
        
    MXFunction([x],[Foo()])
    
    class Foo:
      def __MX__(self):
        return ssym("x")
        
    self.assertRaises(TypeError,lambda : MXFunction([x],[Foo()]))
    
    class Foo:
      def __MX__(self):
        raise Exception("15")
        
    self.assertRaises(TypeError,lambda : MXFunction([x],[Foo()]))

    class Foo:
      pass
        
    self.assertRaises(NotImplementedError,lambda :MXFunction([x],[Foo()]))
    
  def test_OUTPUT(self):
    self.message("OUTPUT typemap")
    a = ssym("A",3,3)
    self.assertTrue(isinstance(qr(a),list))

  def test_cvar(self):
    self.message("We must not have cvar, to avoid bug #652")
    # Wrap all static global things in #ifdef SWIG 
    with self.assertRaises(Exception):
      cvar
      
  def test_ufuncsum(self):
    self.message("ufunc.add")
    
    self.checkarray(DMatrix(sum(DMatrix([1,2,3]))),DMatrix(6))
    
if __name__ == '__main__':
    unittest.main()
