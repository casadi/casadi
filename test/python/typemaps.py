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
import numpy as n
import sys
from numpy import double, int32, ones, matrix, zeros
import unittest
from types import *
from helpers import *
from distutils.version import LooseVersion

scipy_available = True
try:
	from scipy.sparse import *
	from scipy import *
except:
	scipy_available = False
import warnings

print(swig4, systemswig)

class typemaptests(casadiTestCase):

  def setUp(self):
    pass

  def test_memleak(self):
  
   a = np.array([[0, 0]])
   self.assertEqual(sys.getrefcount(a), 2)
   casadi.DM(a)
   self.assertEqual(sys.getrefcount(a), 2)
   casadi.DM(a)
   self.assertEqual(sys.getrefcount(a), 2)
   casadi.SX(a)
   self.assertEqual(sys.getrefcount(a), 2)   
   casadi.SX(a)
   self.assertEqual(sys.getrefcount(a), 2)
   
  def test_0(self):
    self.message("Typemap np.array -> DM")
    arrays = [np.array([[1,2,3],[4,5,6]],dtype=int32),np.array([[1,2,3],[4,5,6]]),np.array([[1,2,3],[4,5,6]],dtype=int), np.array([[1,2],[3,4],[5,6]],dtype=double),np.array([[3.2,4.6,9.9]])]
    for i in range(len(arrays)):
      m = arrays[i]
      zt=c.transpose(c.transpose(m))
      self.assertTrue(isinstance(zt,DM),"DM expected")
      self.checkarray(m,zt,"DM(numpy.ndarray)")
      self.checkarray(m,zt.full(),"DM(numpy.ndarray).full()")
      if scipy_available:
        self.checkarray(m,zt.sparse(),"DM(numpy.ndarray).sparse()")

  def test_1(self):
    self.message("DM -> DM")
    arrays = [DM(Sparsity(4,3,[0,2,2,3],[1,2,1]),[3,2.3,8])]
    for i in range(len(arrays)):
      m = arrays[i]
      zt=c.transpose(c.transpose(m))
      self.assertTrue(isinstance(zt,DM),"DM expected")
      self.checkarray(m,zt,"DM(DM)")
      self.checkarray(m,zt.full(),"DM(DM).full()")
      if scipy_available:
        self.checkarray(m,zt.sparse(),"DM(DM).sparse()")

  def test_2(self):
    self.message("crs_matrix -> DM")
    if not(scipy_available):
      return
    arrays = [csr_matrix( ([3,2.3,8],([0,2,0],[1,1,2])), shape = (3,4), dtype=double ),
              csr_matrix( ([3,2.3,8],([0,2,0],[1,1,2])), shape = (3,4), dtype=int )
              ]
    for i in range(len(arrays)):
      m = arrays[i]
      zt=c.transpose(c.transpose(m))
      self.assertTrue(isinstance(zt,DM),"DM expected")
      self.checkarray(m,zt,"DM(crs_matrix)")
      self.checkarray(m,zt.full(),"DM(crs_matrix).full()")
      if scipy_available:
        self.checkarray(m,zt.sparse(),"DM(crs_matrix).sparse()")


  def test_setget(self):
    return # get/set with return-by-reference has been dropped
    self.message("DM set/get")
    data = n.np.array([3,2.3,8])
    dm=DM(Sparsity(3,4,[0,0,2,3,3],[0,2,0]),[3,2.3,8])

    if scipy_available:
      c=dm.sparse()
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
    self.checkarray(dm.full() > 0,dm,"set(2Dndarray)")
    z=n.matrix(ones((3,4)))
    dm.set(z)
    self.checkarray(dm.full() > 0,dm,"set(2Dmatrix)")
    z=n.zeros((12,5))

    if scipy_available:
      dm.set(c)
      self.checkarray(c,dm,"set(csr_matrix)")

      z=n.zeros(3)
      dm.get_nz(z)
      self.checkarray(n.matrix(z),n.matrix(data),"get(1Dndarray)")
      dm.set_nz(z)

      self.checkarray(c,dm,"set(1Dndarray)")

      dm = dm * 2
      dm.get(c)

      self.checkarray(c,dm,"get(csr_matrix)")

      with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c[0,0]=1
      self.assertRaises(TypeError,lambda :  dm.get(c))

  def test_conversion(self):
    self.message("DM conversions")
    w = DM(Sparsity(4,3,[0,2,2,3],[1,2,1]),[3,2.3,8])
    d = np.array([[1,2,3],[4,5,6]])

    list(w.nonzeros())
    tuple(w.nonzeros())
    w.full()
    np.array(w)
    if check_matrix:
      matrix(w)
    if scipy_available:
      w.sparse()

    self.checkarray(DM(d),d,"DM(numpy.ndarray)")
    #self.checkarray(DM(np.array([1,2,3,4,5,6])),d.ravel(),"DM(numpy.ndarray)")
    #print DM(np.array([1,2,3,4,5,6]))
    #print DM(np.array([1,2,3,6]),2,2).full()

    #print DM(np.array([1,2,3,6]),2,2).full()

  def test_autoconversion(self):
    self.message("Auto conversion DM")
    x=np.array([2.3])
    s = DM([[1,2],[3,4]])
    n = np.array(s)

    self.checkarray(x[0]*s,s*x[0],"")
    self.checkarray(x[0]*s,n*x[0],"")

    self.checkarray(x[0]/s,1/(s/x[0]),"")
    self.checkarray(x[0]/s,x[0]/n,"")

    self.checkarray(x[0]-s,-(s-x[0]),"")
    self.checkarray(x[0]-s,x[0]-n,"")

    self.checkarray(x[0]+s,s+x[0],"")
    self.checkarray(x[0]+s,x[0]+n,"")

    w=np.array([2.3])[0]
    w+=s
    self.checkarray(w,2.3+n,"")

    w=np.array([2.3])[0]
    w-=s
    self.checkarray(w,2.3-n,"")

    w=np.array([2.3])[0]
    w*=s
    self.checkarray(w,2.3*n,"")

    w=np.array([2.3])[0]
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

  @unittest.skipIf(sys.version_info < (3, 5),"Operator was introduced in 3.5")
  def test_matmul(self):
    A = DM([[1,3],[4,5]])
    B = DM([[7,2],[0,9]])
    LA = [A]
    RB = [B]
    if LooseVersion(np.__version__) >= LooseVersion("1.10"):
      LA.append(np.array(A))
      RB.append(np.array(B))
    for L in LA:
      for R in RB:
        #y = L @ R
        y = eval("L @ R",{"L":L,"R":R})
        self.checkarray(y,mtimes(A,B))
    As = SX.sym("x",2,2)
    Bs = SX.sym("x",2,2)
    for L in [As,A]:
      for R in [Bs,B]:
        #y = L @ R
        y = eval("L @ R",{"L":L,"R":R})
        yf = Function('y',[As,Bs],[y])
        self.checkarray(yf(A,B),mtimes(A,B))
    As = MX.sym("x",2,2)
    Bs = MX.sym("x",2,2)
    for L in [As,A]:
      for R in [Bs,B]:
        #y = L @ R
        y = eval("L @ R",{"L":L,"R":R})
        yf = Function('y',[As,Bs],[y])
        self.checkarray(yf(A,B),mtimes(A,B))

  def test_autoconversionMX(self):
    self.message("Auto conversion MX")
    s = DM([[1,2],[3,4]])
    x = SX(3)
    y = MX(3)

    def doit(z,s,fun):
      if type(z) in [type(SX()),type(SX())]:
        ztype = [type(SX()),type(SX())]

      if type(z) in [type(MX())]:
        ztype = [type(MX())]

      r = fun(z,s)

      if type(z) is type(SX()) and type(s) is type(SX()):
        self.assertTrue(type(r) is type(SX()))


      self.assertTrue(type(r) in ztype,"Expected %s but got %s" % (str(ztype),str(type(r))))

      hasNum = True
      if type(s) in [type(SX()),type(MX()),type(SX())]:
        hasNum = False

      if hasNum:
        dummy = [1.3,2.7,9.4,1.0]

        f=Function("f", [z],[r])
        f_in = DM(f.sparsity_in(0),dummy[0:f.nnz_in(0)])
        f_out = f(f_in)

        f_=Function("f", [z],[z])
        f__in = DM(f_.sparsity_in(0),dummy[0:f.nnz_in(0)])
        f__out = f_(f__in)


        self.checkarray(fun(f__out,DM(s)),f_out,"operation",str(f__out)+str(s)+":"+str(fun(f__out,DM(s)))+"->"+str(f_out)+":"+str(s)+str(z)+"->"+str(r))
      else:
        dummy = [1.3,2.7,9.4,1.0]
        dummy2 = [0.3,2.4,1.4,1.7]

        f=Function("f", [z,s],[r])
        f_in = [DM(f.sparsity_in(0),dummy[0:f.nnz_in(0)]), DM(f.sparsity_in(1),dummy2[0:f.nnz_in(1)])]
        f_out = f(*f_in)

        f_=Function("f", [z,s],[z,s])
        f__in = [DM(f_.sparsity_in(0),dummy[0:f.nnz_in(0)]), DM(f_.sparsity_in(1),dummy2[0:f.nnz_in(1)])]
        f__out = f_(*f__in)

        self.checkarray(fun(f__out[0],f__out[1]),f_out,"operation"+str(f__out[0])+","+str(f__out[1])+":"+str(f_out))


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
      doit(z,s,lambda z,s: np.fmin(s,z))
      doit(z,s,lambda z,s: np.fmax(s,z))
      doit(z,s,lambda z,s: np.fmin(s,z))
      doit(z,s,lambda z,s: np.fmax(s,z))
      doit(z,s,lambda z,s: constpow(s,z))
      doit(z,s,lambda z,s: constpow(z,s))
      doit(z,s,lambda z,s: np.arctan2(s,z))
      doit(z,s,lambda z,s: np.arctan2(z,s))
      doit(z,s,lambda z,s: np.copysign(z,s))
      doit(z,s,lambda z,s: np.copysign(s,z))

    nums = [np.array([[1,2],[3,4]]),DM([[1,2],[3,4]]), DM(4), np.array(4),4.0,4]

    ## numeric & SX
    for s in nums:
      for z in [SX.sym("x"), SX.sym("x",2,2)]:
        print("z = %s, s = %s" % (str(z),str(s)))
        print("  z = %s, s = %s" % (type(z),type(s)))
        tests(z,s)

    # numeric & MX
    for s in nums:
      for z in [MX.sym("x",2,2)]:
        print("z = %s, s = %s" % (str(z),str(s)))
        print("  z = %s, s = %s" % (type(z),type(s)))
        tests(z,s)

    ## MX & MX
    for s in [MX.sym("x"),MX.sym("x",2,2)]:
      for z in [MX.sym("x"),MX.sym("x",2,2)]:
        print("z = %s, s = %s" % (str(z),str(s)))
        print("  z = %s, s = %s" % (type(z),type(s)))
        tests(z,s)

    for (s,x,y) in [
                  (np.array([[1,2],[3,4]]),SX.sym("x",2,2),MX.sym("x",2,2))
                  ]:
      for z,ztype in zip([x,y],[[type(SX()),type(SX())],[type(MX())]]):
        print("z = %s, s = %s" % (str(z),str(s)))
        print("  z = %s, s = %s" % (type(z),type(s)))
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
      if type(z) in [type(SX()),type(SX())]:
        ztype = [type(SX()),type(SX())]

      if type(z) in [type(MX())]:
        ztype = [type(MX())]

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

    nums = [np.array([[1,2],[3,4]]),DM([[1,2],[3,4]]), DM(4), np.array(4),4.0,4]

    ## numeric & SX
    for s in nums:
      for z in [SX.sym("x"), SX.sym("x"), SX.sym("x",2,2)]:
        print("z = %s, s = %s" % (str(z),str(s)))
        print("  z = %s, s = %s" % (type(z),type(s)))
        tests(z,s)

    # numeric & MX
    for s in nums:
      for z in [MX.sym("x",2,2)]:
        print("z = %s, s = %s" % (str(z),str(s)))
        print("  z = %s, s = %s" % (type(z),type(s)))
        tests(z,s)

    # SX & SX
    for s in [SX.sym("x"), SX.sym("x"), SX.sym("x",2,2)]:
      for z in [SX.sym("x"),SX.sym("x"), SX.sym("x",2,2)]:
        print("z = %s, s = %s" % (str(z),str(s)))
        print("  z = %s, s = %s" % (type(z),type(s)))
        tests(z,s)

    # MX & MX
    for s in [MX.sym("x"),MX.sym("x",2,2)]:
      for z in [MX.sym("x"),MX.sym("x",2,2)]:
        print("z = %s, s = %s" % (str(z),str(s)))
        print("  z = %s, s = %s" % (type(z),type(s)))
        tests(z,s)

  def test_set(self):
    self.message("DM set on dense matrices")

    # should be integers
    goallist = [1,2,3]
    goal = np.array([goallist]).T

    test={
      "list" : goallist,
      "tuple" : tuple(goallist),
      "array1ddouble" : np.array(goallist,dtype=double),
      "array1dint" : np.array(goallist)
    }
    w=DM(goal)
    self.checkarray(w,goal,"Constructor")

    for name, value in list(test.items()):
      w.nz[:] = value
      self.checkarray(w,goal,"name")

    test={
      "array2ddouble" : np.array([goallist],dtype=double).T,
      "array2dint" : np.array([goallist]).T,
    }
    w=DM(goal)
    self.checkarray(w,goal,"Constructor")

    for name, value in list(test.items()):
      w[:,:] = value
      self.checkarray(w,goal,"name")

  def test_DMSXcast(self):
    self.message("Casting DM to SX")

    W = SX(DM([[1,2,3],[4,5,6]]))

    self.assertEqual(W.size1(),2)
    self.assertEqual(W.size2(),3)

  def test_DMMXcast(self):
    self.message("Casting DM to MX")
    W = MX(DM([[1,2,3],[4,5,6]]))

    self.assertEqual(W.size1(),2)
    self.assertEqual(W.size2(),3)

  def test_DMSX(self):
    self.message("Casting DM to SX")

    w = DM([[1,2,3],[4,5,6]])
    x = SX.sym("x")

    f = Function("f", [x],[w])

    W = f(*f.sx_in())
    self.assertEqual(W.size1(),2)
    self.assertEqual(W.size2(),3)

  def test_DMMX(self):
    self.message("Casting DM to MX")
    w = DM([[1,2,3],[4,5,6]])
    x = MX.sym("x")

    f = Function("f", [x],[w])

    W = f(*f.mx_in())

    self.assertEqual(W.size1(),2)
    self.assertEqual(W.size2(),3)

  def test_setgetslice(self):
    self.message("set/get on DM using slices")

    w = DM([[0,0]])

    A = np.array([[1.0,2],[3,4]])
    B = np.array([[4.0,5],[6,7]])

    w[:,:] = A[0,:]
    self.checkarray(w,A[0,:].reshape((1,-1)),"set")
    B[:1,:] = w
    self.checkarray(B[0,:],A[0,:],"get")

    w = DM([[0],[0]])


    w[:,:] = A[:,0]
    self.checkarray(w,A[:,0],"set")
    B[:,:1] = w
    self.checkarray(B[:,0],A[:,0],"get")

    w = DM([[1,2],[3,4]])
    A = np.zeros((8,7))
    B = np.zeros((8,7))
    A[2:7:3,:7:4] = w
    B[2:7:3,:7:4] = w
    self.checkarray(A,B,"get")

  def test_vertcatprecedence(self):
    self.message("Argument precedence DM")
    a = DM([1,2])
    self.assertTrue(isinstance(vertcat(*[a,a]),DM))

    a = DM([1,2])
    self.assertTrue(isinstance(vertcat(*[a,[1,2,3]]),DM))


    a = MX([1,2])
    print(vertcat(*[a,[1,2,3]]))
    self.assertTrue(isinstance(vertcat(*[a,[1,2,3]]),MX))

  def test_issue190(self):
    self.message("regression test issue #190")
    x=SX.sym("x")
    x * np.array(1)
    x * np.array(1.2)

    SX.sym("x") * np.array(1.0)
    MX.sym("x") * np.array(1.0)

  def test_array_cat(self):
    horzcat(*(SX.sym("x",4,3),np.ones((4,3))))


  def test_issue(self):
    self.message("std::vector<double> typemap.")
    a = np.array([0,2,2,3])
    b = np.array([0.738,0.39,0.99])
    DM(Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.39,0.99])
    DM(Sparsity(4,3,(0,2,2,3),[1,2,1]),[0.738,0.39,0.99])
    DM(Sparsity(4,3,list(a),[1,2,1]),[0.738,0.39,0.99])
    DM(Sparsity(4,3,a,[1,2,1]),[0.738,0.39,0.99])
    DM(Sparsity(4,3,[0,2,2,3],[1,2,1]),(0.738,0.39,0.99))
    DM(Sparsity(4,3,[0,2,2,3],[1,2,1]),list(b))
    DM(Sparsity(4,3,[0,2,2,3],[1,2,1]),b)

  def test_issue314(self):
    self.message("regression test for #314: SX sparsity constructor")
    SX(Sparsity.diag(3),[1,2,3])
  def test_setAll_365(self):
    self.message("ticket #365: DMAtrix.setAll does not work for 1x1 Matrices as input")
    for i in [DM(4),np.array([4]),np.array(4),np.array([[4]]),4,4.0]:
      m = DM.ones(5,5)
      m[:,:] = i
      self.checkarray(m,DM.ones(5,5)*4)

  @unittest.skipIf(sys.version_info>=(3,0), "To lazy to fix")
  def test_issue570(self):
    self.message("Issue #570: long int")
    longint = 10**50
    print(type(longint))
    print(casadi.SX.sym('x') + longint)
    print(longint + casadi.SX.sym('x'))
    print(casadi.SX.sym('x') + longint)
    print(longint + casadi.SX.sym('x'))

  def test_casting_DM(self):
    self.message("casting DM")

    print("cast A")
    x = SX.sym("x")
    f = Function("f", [x],[x])
    class Foo:
      def __DM__(self):
        return DM([4])

    self.assertEqual(f(Foo()),4)
    print("cast B")

    class Foo:
      def __DM__(self):
        return SX([4])


    self.assertRaises(TypeError if systemswig else NotImplementedError,lambda :f(Foo()))
    print("cast C")

    class Foo:
      def __DM__(self):
        raise Exception("15")

    self.assertRaises(TypeError if systemswig else NotImplementedError,lambda :f(Foo()))
    print("cast D")

    class Foo:
      pass

    self.assertRaises(TypeError if systemswig else NotImplementedError,lambda :f(Foo()))
    print("cast E")

  def test_casting_SX(self):
    self.message("casting SX")


    x = SX.sym("x")

    class Foo:
      def __SX__(self):
        return x

    Function("tmp", [x],[Foo()])

    class Foo:
      def __SX__(self):
        return MX.sym("x")

    self.assertRaises(TypeError if systemswig else NotImplementedError,lambda : Function("tmp", [x],[Foo()]))

    class Foo:
      def __SX__(self):
        raise Exception("15")

    self.assertRaises(TypeError if systemswig else NotImplementedError,lambda : Function("tmp", [x],[Foo()]))

    class Foo:
      pass

    self.assertRaises(TypeError if systemswig else NotImplementedError,lambda :Function("tmp", [x],[Foo()]))


  def test_casting_MX(self):
    self.message("casting MX")


    x = MX.sym("x")

    class Foo:
      def __MX__(self):
        return x

    Function("tmp", [x],[Foo()])

    class Foo:
      def __MX__(self):
        return SX.sym("x")

    with self.assertRaises(TypeError if systemswig else NotImplementedError):
      Function("tmp", [x],[Foo()])

    class Foo:
      def __MX__(self):
        raise Exception("15")

    with self.assertRaises(Exception):
      Function("tmp", [x],[Foo()])

    class Foo:
      pass

    with self.assertRaises(TypeError if systemswig else NotImplementedError):
      Function("tmp", [x],[Foo()])

  def test_OUTPUT(self):
    self.message("OUTPUT typemap")
    a = SX.sym("A",3,3)
    
    self.assertTrue(isinstance(qr(a),tuple([tuple]+([list] if swig4 else []))))

  def test_cvar(self):
    self.message("We must not have cvar, to avoid bug #652")
    # Wrap all static global things in #ifdef SWIG
    with self.assertRaises(Exception):
      cvar

  def test_ufuncsum(self):
    self.message("ufunc.add")

    self.checkarray(DM(np.sum(DM([1,2,3]).nonzeros())),DM(6))

  def test_sxmatrix(self):

    def val(a):
      f = Function("f", [],[a])
      f_out = f.call([])
      return f_out[0]

    for i in [SX(1),1,1.0]:
      a = np.array([[1,2],[3,4]])
      print(val(SX(a)))
      print(val(SX(a.T)))

      self.checkarray(val(SX(a)),DM([[1,2],[3,4]]))
      self.checkarray(val(SX(a.T).T),DM([[1,2],[3,4]]))

      if check_matrix:
        a = numpy.matrix([[1,2],[3,4]])

        print(val(SX(a)))
        print(DM([[1,2],[3,4]]))

        self.checkarray(val(SX(a)),DM([[1,2],[3,4]]))
        self.checkarray(val(SX(a.T).T),DM([[1,2],[3,4]]))

  def test_issue1158(self):
    A = numpy.zeros((0,2))
    a = DM(A)
    self.assertEqual(a.shape[0],0)
    self.assertEqual(a.shape[1],2)

  def test_matrices(self):

    Ds = ([
          numpy.matrix([[1,2],[3,4]]),
          numpy.matrix([[1,2],[3,4.0]])] if check_matrix else [])+[
          np.array([[1,2],[3,4]]),
          np.array([[1,2],[3,4.0]]),
        ]

    if scipy_available:
      Ds+=[
          csc_matrix(([1.0,3.0,2.0,4.0],[0,1,0,1],[0,2,4]),shape=(2,2),dtype=numpy.double),
          csc_matrix(([1,3,2,4],[0,1,0,1],[0,2,4]),shape=(2,2),dtype=numpy.int),
          DM([[1,2],[3,4]]).sparse()
      ]


    for D in Ds:
      print(D)
      d = DM.ones(2,2)

      x = SX.sym("x",d.sparsity())
      f = Function("f", [x],[x])
      fin = DM(d.sparsity(),0)
      fin[:,:] = D

      self.checkarray(fin,DM([[1,2],[3,4]]))
      d[:,:] = D
      self.checkarray(d,DM([[1,2],[3,4]]))

  def test_issue1217(self):
    a = np.array([0,SX.sym("x")])

    print(if_else(0,a,a))

  def test_issue1373(self):
    print(np.array(casadi.DM([2])))
    print(np.array(casadi.DM([1,2,3.0])))

  def test_None(self):
    #self.assertFalse(None==DM(3))
    b = np.atleast_2d(None)
    with self.assertRaises(TypeError if systemswig else NotImplementedError):
      c = repmat(b, 1, 1)

  @requires_nlpsol("ipopt")
  def testGenericTypeBoolean(self):
    x=SX.sym("x")
    with self.assertRaises(RuntimeError):
      solver = nlpsol("mysolver", "ipopt", {"x":x,"f":x**2}, {"ipopt": {"acceptable_tol": SX.sym("x")}})

    nlpsol("mysolver", "ipopt", {"x":x,"f":x**2}, {"ipopt": {"acceptable_tol": 1}})
  
  """
  def test_to_longlong(self):
    a = IM(10)


    b = a**15

    self.assertEqual(int(b),10**15)
  """

  def test_buglonglong(self):
    x = SX.sym("x")

    jacobian(x/1.458151064450277e-12,x)


  def test_issue_2625(self):
    # This is obviously a bug
    self.checkarray(np.inner(DM([1,0,1]),DM([1,0,1])), np.array([[1,0,1],[0,0,0],[1,0,1]]))
    self.checkarray(np.outer(DM([1,0,1]),DM([1,0,1])), np.array([[1,0,1],[0,0,0],[1,0,1]]))
    #print(np.logical_or(DM([1,0,1]),DM([1,0,12])))
    self.checkarray(np.add.accumulate(DM([1,0,1])),np.array([[1],[1],[2]]))
    self.checkarray(np.cumsum(DM([1,0,1])), np.array([1,1,2]))
    self.checkarray(np.sum(DM([1,0,1])),np.array(2))
    self.assertFalse(np.all(DM([1,0,1])))
    self.assertTrue(np.any(DM([1,0,1])))
    self.checkarray(np.sum(DM([[1,0,1],[0,1,0]])),np.array(3))
    self.assertFalse(np.all(DM([[1,0,1],[0,1,0]])))
    self.assertTrue(np.any(DM([[1,0,1],[0,1,0]])))

    for M in [SX,MX]:
      with self.assertRaises(Exception):
        np.inner(M([1,0,1]),M([1,0,1]))
      with self.assertRaises(Exception):
        np.outer(M([1,0,1]),M([1,0,1]))

      with self.assertRaises(Exception):
        np.add.accumulate(M([1,0,1]))

      with self.assertRaises(Exception):
        np.cumsum(M([1,0,1]))

      with self.assertRaises(Exception):
        np.sum(M([1,0,1]))
      with self.assertRaises(Exception):
        np.all(M([1,0,1]))
      with self.assertRaises(Exception):
        np.any(M([1,0,1]))
      with self.assertRaises(Exception):
        np.sum(M([[1,0,1],[0,1,0]]))
      with self.assertRaises(Exception):
        np.all(M([[1,0,1],[0,1,0]]))
      with self.assertRaises(Exception):
        np.any(M([[1,0,1],[0,1,0]]))

if __name__ == '__main__':
    unittest.main()
