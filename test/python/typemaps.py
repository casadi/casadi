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
import casadi as c
import numpy
import numpy as n
import numpy as np
import sys
from numpy import double, int32, ones, matrix, zeros
import unittest
import math
from types import *
from helpers import *
from looseversion import LooseVersion

scipy_available = True
try:
	from scipy.sparse import *
	from scipy import *
except:
	scipy_available = False
import warnings

print(swig4, systemswig)

class typemaptests(casadiTestCase):

  def _assert_silent_numeric(self, val):
    # casadi.<fn>(numpy) is the inverse direction: it never consults numpy_mode
    # (only the numpy.foo(casadi) gateway does), so the SWIG typemap densifies
    # the ndarray to a casadi DM, the same in every mode and every numpy version.
    self.assertIsInstance(val, ca.DM)

  def setUp(self):
    pass

  def test_memleak(self):
  
   a = np.array([[0, 0]])
   r = sys.getrefcount(a)
   ca.DM(a)
   self.assertEqual(sys.getrefcount(a), r)
   ca.DM(a)
   self.assertEqual(sys.getrefcount(a), r)
   ca.SX(a)
   self.assertEqual(sys.getrefcount(a), r)   
   ca.SX(a)
   self.assertEqual(sys.getrefcount(a), r)
   
  def test_sanity(self):
    a = ca.DM([1,2,3])
    if abs(ca.norm_inf(a))>1:
        self.assertTrue(True)
    else:
        self.assertTrue(False)
   
    a = ca.SX(1)
    if a>0:
        self.assertTrue(True)
    else:
        self.assertTrue(False)
    if abs(a)>0:
        self.assertTrue(True)
    else:
        self.assertTrue(False)
    a = ca.MX(1)
    if a>0:
        self.assertTrue(True)
    else:
        self.assertTrue(False)
    if abs(a)>0:
        self.assertTrue(True)
    else:
        self.assertTrue(False)
   
  def test_0(self):
    self.message("Typemap np.array -> DM")
    arrays = [np.array([[1,2,3],[4,5,6]],dtype=int32),np.array([[1,2,3],[4,5,6]]),np.array([[1,2,3],[4,5,6]],dtype=int), np.array([[1,2],[3,4],[5,6]],dtype=double),np.array([[3.2,4.6,9.9]])]
    for i in range(len(arrays)):
      m = arrays[i]
      zt=c.transpose(c.transpose(m))
      self.assertTrue(isinstance(zt,ca.DM),"DM expected")
      self.checkarray(m,zt,"DM(numpy.ndarray)")
      self.checkarray(m,zt.full(),"DM(numpy.ndarray).full()")
      if scipy_available:
        self.checkarray(m,zt.sparse(),"DM(numpy.ndarray).sparse()")

  def test_1(self):
    self.message("DM -> DM")
    arrays = [ca.DM(ca.Sparsity(4,3,[0,2,2,3],[1,2,1]),[3,2.3,8])]
    for i in range(len(arrays)):
      m = arrays[i]
      zt=c.transpose(c.transpose(m))
      self.assertTrue(isinstance(zt,ca.DM),"DM expected")
      self.checkarray(m,zt,"DM(DM)")
      self.checkarray(m,zt.full(),"DM(DM).full()")
      if scipy_available:
        self.checkarray(m,zt.sparse(),"DM(DM).sparse()")

  def test_2(self):
    self.message("crs_matrix -> DM")
    if not(scipy_available):
      return
    arrays = [csr_matrix( ([3,2.3,8],([0,2,0],[1,1,2])), shape = (3,4), dtype=double ),  # pyright: ignore[reportUnboundVariable]
              csr_matrix( ([3,2.3,8],([0,2,0],[1,1,2])), shape = (3,4), dtype=int )  # pyright: ignore[reportUnboundVariable]
              ]
    for i in range(len(arrays)):
      m = arrays[i]
      zt=c.transpose(c.transpose(m))
      self.assertTrue(isinstance(zt,ca.DM),"DM expected")
      self.checkarray(m,zt,"DM(crs_matrix)")
      self.checkarray(m,zt.full(),"DM(crs_matrix).full()")
      if scipy_available:
        self.checkarray(m,zt.sparse(),"DM(crs_matrix).sparse()")


  def test_setget(self):
    return # get/set with return-by-reference has been dropped
    self.message("DM set/get")
    data = n.np.array([3,2.3,8])
    dm=ca.DM(ca.Sparsity(3,4,[0,0,2,3,3],[0,2,0]),[3,2.3,8])

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
    w = ca.DM(ca.Sparsity(4,3,[0,2,2,3],[1,2,1]),[3,2.3,8])
    d = np.array([[1,2,3],[4,5,6]])

    list(w.nonzeros())
    tuple(w.nonzeros())
    w.full()
    np.array(w)
    if check_matrix:
      matrix(w)
    if scipy_available:
      w.sparse()

    self.checkarray(ca.DM(d),d,"DM(numpy.ndarray)")
    #self.checkarray(DM(np.array([1,2,3,4,5,6])),d.ravel(),"DM(numpy.ndarray)")
    #print DM(np.array([1,2,3,4,5,6]))
    #print DM(np.array([1,2,3,6]),2,2).full()

    #print DM(np.array([1,2,3,6]),2,2).full()

  def test_autoconversion(self):
    self.message("Auto conversion DM")
    x=np.array([2.3])
    s = ca.DM([[1,2],[3,4]])
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
    A = ca.DM([[1,3],[4,5]])
    B = ca.DM([[7,2],[0,9]])
    LA = [A]  # type: list
    RB = [B]  # type: list
    if LooseVersion(np.__version__) >= LooseVersion("1.10"):
      LA.append(np.array(A))
      RB.append(np.array(B))
    for L in LA:
      for R in RB:
        #y = L @ R
        y = eval("L @ R",{"L":L,"R":R})
        self.checkarray(y,A @ B)
    As = ca.SX.sym("x",2,2)
    Bs = ca.SX.sym("x",2,2)
    for L in [As,A]:
      for R in [Bs,B]:
        #y = L @ R
        y = eval("L @ R",{"L":L,"R":R})
        yf = ca.Function('y',[As,Bs],[y])
        self.checkarray(yf(A,B),A @ B)
    As = ca.MX.sym("x",2,2)
    Bs = ca.MX.sym("x",2,2)
    for L in [As,A]:
      for R in [Bs,B]:
        #y = L @ R
        y = eval("L @ R",{"L":L,"R":R})
        yf = ca.Function('y',[As,Bs],[y])
        self.checkarray(yf(A,B),A @ B)

  def test_autoconversionMX(self):
    self.message("Auto conversion MX")
    s = ca.DM([[1,2],[3,4]])
    x = ca.SX(3)
    y = ca.MX(3)

    def doit(z,s,fun):
      if type(z) in [type(ca.SX()),type(ca.SX())]:
        ztype = [type(ca.SX()),type(ca.SX())]

      if type(z) in [type(ca.MX())]:
        ztype = [type(ca.MX())]

      r = fun(z,s)

      if type(z) is type(ca.SX()) and type(s) is type(ca.SX()):
        self.assertTrue(type(r) is type(ca.SX()))


      self.assertTrue(type(r) in ztype,"Expected %s but got %s" % (str(ztype),str(type(r))))

      hasNum = True
      if type(s) in [type(ca.SX()),type(ca.MX()),type(ca.SX())]:
        hasNum = False

      if hasNum:
        dummy = [1.3,2.7,9.4,1.0]

        f=ca.Function("f", [z],[r])
        f_in = ca.DM(f.sparsity_in(0),dummy[0:f.nnz_in(0)])
        f_out = f(f_in)

        f_=ca.Function("f", [z],[z])
        f__in = ca.DM(f_.sparsity_in(0),dummy[0:f.nnz_in(0)])
        f__out = f_(f__in)


        self.checkarray(fun(f__out,ca.DM(s)),f_out,"operation",str(f__out)+str(s)+":"+str(fun(f__out,ca.DM(s)))+"->"+str(f_out)+":"+str(s)+str(z)+"->"+str(r))
      else:
        dummy = [1.3,2.7,9.4,1.0]
        dummy2 = [0.3,2.4,1.4,1.7]

        f=ca.Function("f", [z,s],[r])
        f_in = [ca.DM(f.sparsity_in(0),dummy[0:f.nnz_in(0)]), ca.DM(f.sparsity_in(1),dummy2[0:f.nnz_in(1)])]
        f_out = f(*f_in)

        f_=ca.Function("f", [z,s],[z,s])
        f__in = [ca.DM(f_.sparsity_in(0),dummy[0:f.nnz_in(0)]), ca.DM(f_.sparsity_in(1),dummy2[0:f.nnz_in(1)])]
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
      # casadi fmin/fmax keep autoconversion casadi-typed; the numpy
      # spellings np.fmin/np.fmax now route through the opt-in casadi-aware
      # numpy support (issue #2959) and are covered in numpy_interop.py.
      doit(z,s,lambda z,s: ca.fmin(s,z))
      doit(z,s,lambda z,s: ca.fmax(s,z))
      doit(z,s,lambda z,s: ca.fmin(s,z))
      doit(z,s,lambda z,s: ca.fmax(s,z))
      doit(z,s,lambda z,s: ca.constpow(s,z))
      doit(z,s,lambda z,s: ca.constpow(z,s))
      doit(z,s,lambda z,s: ca.arctan2(s,z))
      doit(z,s,lambda z,s: ca.arctan2(z,s))
      doit(z,s,lambda z,s: ca.copysign(z,s))
      doit(z,s,lambda z,s: ca.copysign(s,z))
      doit(z,s,lambda z,s: z @ s)
      doit(z,s,lambda z,s: s @ z)
      # Py3 numeric dunders added alongside the operator audit: in-place
      # arith, divmod, pow-with-modulo, and the math protocol unary ops.
      doit(z,s,lambda z,s: z.__iadd__(s))
      doit(z,s,lambda z,s: z.__isub__(s))
      doit(z,s,lambda z,s: z.__imul__(s))
      doit(z,s,lambda z,s: z.__itruediv__(s))
      doit(z,s,lambda z,s: z.__ifloordiv__(s))
      doit(z,s,lambda z,s: z.__imod__(s))
      doit(z,s,lambda z,s: z.__ipow__(s))
      doit(z,s,lambda z,s: divmod(z,s)[0])
      doit(z,s,lambda z,s: divmod(z,s)[1])
      doit(z,s,lambda z,s: divmod(s,z)[0])
      doit(z,s,lambda z,s: divmod(s,z)[1])
      doit(z,s,lambda z,s: pow(z, s, 5.0))
      doit(z,s,lambda z,s: round(z))
      doit(z,s,lambda z,s: math.trunc(z))
      doit(z,s,lambda z,s: math.floor(z))
      doit(z,s,lambda z,s: math.ceil(z))

    nums = [np.array([[1,2],[3,4]]),ca.DM([[1,2],[3,4]]), ca.DM(4), np.array(4),4.0,4]

    ## numeric & SX
    for s in nums:
      for z in [ca.SX.sym("x"), ca.SX.sym("x",2,2)]:
        print("z = %s, s = %s" % (str(z),str(s)))
        print("  z = %s, s = %s" % (type(z),type(s)))
        tests(z,s)

    # numeric & MX
    for s in nums:
      for z in [ca.MX.sym("x",2,2)]:
        print("z = %s, s = %s" % (str(z),str(s)))
        print("  z = %s, s = %s" % (type(z),type(s)))
        tests(z,s)

    # SX & SX
    for s in [ca.SX.sym("x"), ca.SX.sym("x"), ca.SX.sym("x",2,2)]:
      for z in [ca.SX.sym("x"),ca.SX.sym("x"), ca.SX.sym("x",2,2)]:
        print("z = %s, s = %s" % (str(z),str(s)))
        print("  z = %s, s = %s" % (type(z),type(s)))
        tests(z,s)

    ## MX & MX
    for s in [ca.MX.sym("x"),ca.MX.sym("x",2,2)]:
      for z in [ca.MX.sym("x"),ca.MX.sym("x",2,2)]:
        print("z = %s, s = %s" % (str(z),str(s)))
        print("  z = %s, s = %s" % (type(z),type(s)))
        tests(z,s)

    for (s,x,y) in [
                  (np.array([[1,2],[3,4]]),ca.SX.sym("x",2,2),ca.MX.sym("x",2,2))
                  ]:
      for z,ztype in zip([x,y],[[type(ca.SX()),type(ca.SX())],[type(ca.MX())]]):
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
        
  def test_numpy_preserve_type(self):
    """issue #2959: GlobalOptions.setNumpyMode is the sole gateway to the
    casadi-aware numpy support.  mode 0 (default): an explicit numpy.foo(M)
    on a casadi value uses legacy 3.7.2 behaviour + a FutureWarning (numeric
    densifies to a numpy result; symbolic returns a casadi value).  mode -1:
    the same, silently.  mode 1: numpy.foo(M) returns a casadi.ArrayInterface
    following numpy's shape/axis contract, for DM/SX/MX alike.  Operator
    arithmetic (M + x, numpy_scalar - M) and casadi.foo(numpy) always keep
    their casadi-typed 3.7.2 behaviour, in every mode."""
    Arr = c.ArrayInterface
    self.assertTrue(hasattr(ca.GlobalOptions, "setNumpyMode"))
    self.assertEqual(ca.GlobalOptions.getNumpyMode(), 0)
    try:
      M = ca.DM([[1.0, 2.0], [3.0, 4.0]])

      # --- default: explicit numpy.foo densifies + warns (numeric).
      # simplefilter("always") overrides the once-per-process default filter,
      # so every call warns and each op is checked. ---
      for label, fn in [("np.sin", lambda x: numpy.sin(x)),
                        ("np.sum", lambda x: numpy.sum(x)),
                        ("np.reshape", lambda x: numpy.reshape(x, (1, 4)))]:
        with warnings.catch_warnings(record=True) as w:
          warnings.simplefilter("always")
          r = fn(M)
          self.assertFalse(isinstance(r, (ca.DM, ca.SX, ca.MX)),
                          "%s default must densify, got %s" %
                          (label, type(r).__name__))
          self.assertTrue(any(issubclass(x.category, FutureWarning)
                              for x in w), "%s must warn" % label)
          self.assertTrue(any("setNumpyMode" in str(x.message)
                              for x in w), "%s must point at the flag" % label)

      # --- operator arithmetic stays casadi-typed and silent ---
      with warnings.catch_warnings():
        warnings.simplefilter("error", category=FutureWarning)
        self.assertTrue(isinstance(numpy.float64(2.0) - M, ca.DM))
        self.assertTrue(isinstance(M * numpy.array([[1.0, 1.0], [1.0, 1.0]]), ca.DM))
        for cls in [ca.SX, ca.MX]:
          s = cls.sym("s", 2, 2)
          self.assertTrue(isinstance(numpy.float64(2.0) - s, cls))

      # --- default + symbolic explicit numpy.foo: 3.7.2 returned a casadi
      #     value, so warn + return the casadi type (np.sin(MX) -> MX) ---
      for cls in [ca.SX, ca.MX]:
        s = cls.sym("s", 2, 2)
        with warnings.catch_warnings(record=True) as w:
          warnings.simplefilter("always")
          r = numpy.sin(s)
          self.assertIsInstance(r, cls)
          self.assertTrue(any(issubclass(x.category, FutureWarning)
                              for x in w), "np.sin(%s) must warn" % cls.__name__)

      # --- casadi.foo(numpy) keeps exact 3.7.2 behaviour in EVERY mode:
      #     it densifies the numpy array to a casadi DM silently (no warning).
      #     This is the inverse direction; it never consults numpy_mode, so it
      #     is identical across modes.  casadi values / scalars are unaffected. ---
      arr = numpy.array([1.0, 2.0, 3.0])
      with warnings.catch_warnings():
        warnings.simplefilter("error", category=FutureWarning)
        self._assert_silent_numeric(ca.sin(arr))                    # numpy -> DM, silent
        self.assertTrue(isinstance(ca.sin(ca.DM([1.0, 2.0])), ca.DM))   # casadi: silent
        self.assertTrue(isinstance(ca.sin(ca.SX.sym("z")), ca.SX))      # casadi: silent
        self.assertEqual(type(ca.sin(0.3)), float)                # scalar: silent

      # --- mode -1: legacy behaviour, but SILENT (no FutureWarning) ---
      ca.GlobalOptions.setNumpyMode(-1)
      with warnings.catch_warnings():
        warnings.simplefilter("error", category=FutureWarning)
        r = numpy.sin(M)                       # densifies, no warning
        self.assertFalse(isinstance(r, (ca.DM, ca.SX, ca.MX)))
        self.assertIsInstance(numpy.sin(ca.SX.sym("z", 2, 2)), ca.SX)  # symbolic stays casadi
        self._assert_silent_numeric(ca.sin(numpy.array([1.0, 2.0])))  # inverse: silent
      ca.GlobalOptions.setNumpyMode(0)

      # --- preserve: returns ArrayInterface, numpy shape/axis contract ---
      ca.GlobalOptions.setNumpyMode(1)
      with warnings.catch_warnings():
        warnings.simplefilter("error", category=FutureWarning)
        r = numpy.sin(M)
        self.assertIsInstance(r, Arr)
        self.assertEqual(r.shape, (2, 2))
        tot = numpy.sum(M)                  # full reduction -> scalar
        self.assertIsInstance(tot, Arr)
        self.assertEqual(float(tot), 10.0)
        col = numpy.sum(M, axis=0)          # axis 0 dropped -> (2,)
        self.assertEqual(col.shape, (2,))
        # symbolic flows through the same gateway
        sx = ca.SX.sym("s", 2, 2)
        e = numpy.sin(sx)
        self.assertIsInstance(e, Arr)
        # numpy.sin's stub return type is NDArray (it can't see our
        # __array_ufunc__), so pyright doesn't know `e` is an ArrayInterface.
        f = ca.Function("f", [sx], [e.to_casadi()])  # pyright: ignore[reportAttributeAccessIssue]
        self.checkarray(f(M), ca.DM(numpy.sin(numpy.array(M))), "symbolic np.sin")
        # casadi.sin(numpy) is unchanged even under mode 1: silent
        self._assert_silent_numeric(ca.sin(numpy.array([1.0, 2.0, 3.0])))
    finally:
      ca.GlobalOptions.setNumpyMode(0)

  def test_numpy_mode_behaviour(self):
    """issue #2959 behavioural reference -- doubles as release notes.

    A PLAIN NUMPY ARRAY argument behaves exactly like casadi 3.7.2 in EVERY
    mode; no FutureWarning is ever raised:

        a = np.ones((2, 3))
        ca.sum(a) -> DM             np.sum(a) -> np.float64
        ca.sin(a) -> DM             np.sin(a) -> np.ndarray

    The behaviour diverges ONLY for a casadi value, e.g.
    x = ca.MX.sym("x", ca.Sparsity.upper(3)):

        ca.sum(x), ca.sin(x) -> MX always (casadi-native, silent)
        np.sum(x), np.sin(x) -> mode -1: MX, silent
                                mode  0: MX + FutureWarning   (default)
                                mode  1: ca.ArrayInterface
    """
    import numpy as np
    Arr = c.ArrayInterface
    orig = ca.GlobalOptions.getNumpyMode()
    try:
      # (1) numpy-array argument: identical in every mode, and never warns.
      a = np.ones((2, 3))
      for mode in (-1, 0, 1):
        ca.GlobalOptions.setNumpyMode(mode)
        with warnings.catch_warnings():
          warnings.simplefilter("error", category=FutureWarning)
          self.assertIsInstance(c.sum(a), ca.DM)                 # casadi densifies -> DM
          self.assertEqual(float(c.sum(a)), 6.0)
          self.assertIsInstance(np.sum(a), np.floating)       # pure numpy
          self.assertEqual(float(np.sum(a)), 6.0)
          self.assertIsInstance(c.sin(a), ca.DM)                 # casadi densifies -> DM
          self.assertIsInstance(np.sin(a), np.ndarray)        # pure numpy

      # (2) casadi-value argument: ca.* stays casadi-native and silent always.
      x = ca.MX.sym("x", ca.Sparsity.upper(3))
      for mode in (-1, 0, 1):
        ca.GlobalOptions.setNumpyMode(mode)
        with warnings.catch_warnings():
          warnings.simplefilter("error", category=FutureWarning)
          self.assertIsInstance(c.sum(x), ca.MX)
          self.assertIsInstance(c.sin(x), ca.MX)

      # (3) numpy.* on a casadi value: the one place the mode matters.
      ca.GlobalOptions.setNumpyMode(-1)                          # legacy, silent
      with warnings.catch_warnings():
        warnings.simplefilter("error", category=FutureWarning)
        self.assertIsInstance(np.sum(x), ca.MX)
        self.assertIsInstance(np.sin(x), ca.MX)

      ca.GlobalOptions.setNumpyMode(0)                           # legacy + FutureWarning
      with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        self.assertIsInstance(np.sum(x), ca.MX)
        self.assertTrue(any(issubclass(e.category, FutureWarning) for e in w),
                        "np.sum(casadi) must warn in mode 0")
      with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        self.assertIsInstance(np.sin(x), ca.MX)
        self.assertTrue(any(issubclass(e.category, FutureWarning) for e in w),
                        "np.sin(casadi) must warn in mode 0")
      # (the "fires once by default, -W always shows all" property is checked
      # in test_numpy_warning_shown_once)

      ca.GlobalOptions.setNumpyMode(1)                           # numpy-aware
      with warnings.catch_warnings():
        warnings.simplefilter("error", category=FutureWarning)
        self.assertIsInstance(np.sum(x), Arr)
        self.assertIsInstance(np.sin(x), Arr)
    finally:
      ca.GlobalOptions.setNumpyMode(orig)

  def test_numpy_warning_shown_once(self):
    """The legacy-mode FutureWarning is deduped to once per process by a
    lowest-precedence 'once' filter -- so a user is never flooded WITHOUT
    touching the warnings module.  An AI-aided developer surfaces every
    occurrence with the STANDARD `python -W always` (or PYTHONWARNINGS=
    always), and `-W error` makes it fatal -- no casadi-specific knob."""
    import subprocess, sys, os
    code = ("import casadi as ca, numpy as np\n"
            "M = ca.DM([[1.,2.],[3.,4.]])\n"
            "x = ca.MX.sym('x', 2, 2)\n"
            "for _ in range(20):\n"
            "    np.sum(M); np.sin(x)\n")
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(sys.path))
    def run(*flags):
      return subprocess.run([sys.executable, *flags, "-c", code],
                            capture_output=True, text=True, env=env)
    # count the unique notice text (the echoed source line also contains the
    # word "FutureWarning", so don't count that).
    needle = "casadi: a numpy function"
    default = run()
    self.assertEqual(default.stderr.count(needle), 1,
                     "default filters: exactly one notice\n" + default.stderr)
    every = run("-W", "always")
    self.assertGreater(every.stderr.count(needle), 1,
                       "-W always: every occurrence\n" + every.stderr)
    fatal = run("-W", "error::FutureWarning")
    self.assertNotEqual(fatal.returncode, 0,
                        "-W error::FutureWarning must make it fatal")

  def test_numpy_no_global_filter(self):
    """issue #2959: the once-per-process dedup uses the stdlib default action
    (a constant notice emitted from one fixed line), NOT a global warnings
    filter -- so importing casadi must not mutate warnings.filters for our
    notice, leaving other packages' warnings untouched."""
    mine = [f for f in warnings.filters
            if f[1] is not None and "numpy function" in str(f[1])]
    self.assertEqual(mine, [], "issue #2959 must not install a warnings filter")

  def test_numpy_unary_ufuncs_warn(self):
    """issue #2959: np.abs / np.negative / np.positive are explicit np.foo,
    so they warn and follow the type rule (DM -> numpy, SX/MX -> casadi).
    The Python operators abs(M) / -M and the binary-operator ufuncs (np.add)
    stay silent casadi, because those are reachable via operator interop."""
    import numpy as np
    orig = ca.GlobalOptions.getNumpyMode()
    ca.GlobalOptions.setNumpyMode(0)
    try:
      D = ca.DM([[1.0, -2.0]]); S = ca.SX.sym("s", 1, 2); M = ca.MX.sym("m", 1, 2)
      for fn in (np.abs, np.negative, np.positive):
        with warnings.catch_warnings(record=True) as w:
          warnings.simplefilter("always")
          rd = fn(D)
          self.assertFalse(isinstance(rd, (ca.DM, ca.SX, ca.MX)),  # numeric -> numpy
                           "%s(DM) must densify" % fn.__name__)
          self.assertTrue(any(issubclass(e.category, FutureWarning) for e in w),
                          "%s(DM) must warn" % fn.__name__)
        for sym in (S, M):
          with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            r = fn(sym)
            self.assertIsInstance(r, type(sym))           # symbolic -> casadi
            self.assertTrue(any(issubclass(e.category, FutureWarning) for e in w),
                            "%s(%s) must warn" % (fn.__name__, type(sym).__name__))
      # operators and binary-operator ufuncs stay silent casadi (no warning)
      with warnings.catch_warnings():
        warnings.simplefilter("error", category=FutureWarning)
        self.assertIsInstance(abs(D), ca.DM)         # builtin abs -> __abs__
        self.assertIsInstance(-D, ca.DM)             # __neg__
        self.assertIsInstance(np.add(D, D), ca.DM)   # operator ufunc
    finally:
      ca.GlobalOptions.setNumpyMode(orig)

  def test_operator_array_priority_deferral(self):
    """issue #2959 (option A): a casadi binary operator yields to an operand
    with a higher __array_priority__, so an ArrayInterface wins in BOTH
    directions -- `w op M` and `M op w` both return an ArrayInterface of the
    backing type -- while ordinary casadi / numpy / scalar interop is
    unchanged."""
    import numpy as np
    Arr = c.ArrayInterface
    ops = [lambda a, b: a + b,  lambda a, b: a - b,  lambda a, b: a * b,
           lambda a, b: a / b,  lambda a, b: a // b, lambda a, b: a % b,
           lambda a, b: a @ b,  lambda a, b: a ** b, lambda a, b: a < b,
           lambda a, b: a <= b, lambda a, b: a > b,  lambda a, b: a >= b,
           lambda a, b: a == b, lambda a, b: a != b]
    for name, raw in [("DM", ca.DM([[1.0, 2.0], [3.0, 4.0]])),
                      ("SX", ca.SX.sym("s", 2, 2)),
                      ("MX", ca.MX.sym("m", 2, 2))]:
      w = Arr(raw)
      expected = "ArrayInterface" + name
      for op in ops:
        self.assertEqual(type(op(w, raw)).__name__, expected)   # wrapper on left
        self.assertEqual(type(op(raw, w)).__name__, expected)   # raw on left -> defers
    # interop unchanged: casadi promotion, casadi-beats-ndarray, scalars
    D = ca.DM([[1.0, 2.0]]); M = ca.MX.sym("m", 1, 2); a = np.array([[3.0, 4.0]])
    self.assertIsInstance(D + M, ca.MX)         # DM/MX promotion
    self.assertIsInstance(M + D, ca.MX)
    self.assertIsInstance(D + a, ca.DM)         # DM on left beats ndarray
    self.assertIsInstance(D + 2, ca.DM)         # python scalar
    self.assertIsInstance(2 * D, ca.DM)

  def test_array_alias(self):
    """issue #2959: casadi.array is a qualified-only lowercase factory for the
    array-semantics wrapper (like pd.array / jnp.array).  It must NOT leak into
    `from casadi import *` -- a bare `array` would shadow numpy.array."""
    import numpy as np
    self.assertIs(ca.array, ca.ArrayInterface)
    self.assertIsInstance(ca.array(ca.SX.sym("s", 2, 2)), ca.ArrayInterface)
    # star-import must NOT bring `array` ...
    ns = {}
    exec("from casadi import *", ns)
    self.assertNotIn("array", ns)
    self.assertIn("ArrayInterface", ns)         # ... while ArrayInterface still does
    # ... so numpy.array survives the common import combo
    ns2 = {}
    exec("from numpy import *\nfrom casadi import *", ns2)
    self.assertIs(ns2["array"], np.array)

  def test_squeeze(self):
    """issue #2959: numpy .squeeze (drop length-1 axes) on the ArrayInterface."""
    import numpy as np
    Arr = ca.ArrayInterface
    for shp in [(3, 1), (1, 3), (1, 1), (2, 3), (1, 4), (4, 1)]:
      A = np.arange(float(np.prod(shp))).reshape(shp)
      s = Arr(ca.DM(A)).squeeze()
      self.assertEqual(s.shape, A.squeeze().shape)
      self.checkarray(np.asarray(s).reshape(A.squeeze().shape), A.squeeze(), "squeeze %s" % (shp,))
    self.assertEqual(Arr(ca.DM(np.ones((1, 3)))).squeeze(0).shape, (3,))   # axis arg
    with self.assertRaises(ValueError):
      Arr(ca.DM(np.ones((1, 3)))).squeeze(1)                                # non-1 axis
    self.assertIsInstance(Arr(ca.SX.sym("s", 1, 3)).squeeze(), Arr)        # symbolic

  def test_flat_flatten_ravel(self):
    """issue #2959: numpy .flat iterator (C-order iterate / index / slice /
    write-back) plus .flatten() and .ravel() on the ArrayInterface."""
    import numpy as np
    Arr = ca.ArrayInterface
    A = np.arange(6.0).reshape(2, 3)        # C-order: 0..5
    x = Arr(ca.DM(A))
    self.assertEqual(len(x.flat), A.size)
    self.assertEqual([float(v) for v in x.flat], list(A.flat))        # iterate
    self.assertEqual(float(x.flat[4]), A.flat[4])                     # int index
    self.assertEqual([float(v) for v in x.flat[1:5]], list(A.flat[1:5]))  # slice
    self.assertEqual([float(v) for v in np.asarray(x.flatten())], list(A.flatten()))
    self.assertEqual([float(v) for v in np.asarray(x.ravel())], list(A.ravel()))
    # write-back through .flat mutates the owner, in C-order
    x.flat[2] = 99.0;  A.flat[2] = 99.0
    x.flat[0:2] = -1;  A.flat[0:2] = -1
    self.checkarray(np.asarray(x.__array__()), A, "flat setitem")
    # symbolic: .flat yields wrappers, .flatten stays casadi-typed
    s = Arr(ca.SX.sym("s", 2, 2))
    self.assertEqual([type(e).__name__ for e in s.flat], ["ArrayInterfaceSX"] * 4)
    self.assertIsInstance(s.flatten(), Arr)

  def test_experimental_numpy_array(self):
    """issue #2959: the python-only numpy-semantics array wrapper
    (casadi.ArrayInterface).  Variable ndim (0/1/2), numpy row indexing /
    iteration, right-aligned broadcasting, matmul, and numpy dispatch that
    returns wrapper instances -- all on top of the 2-D casadi values, for
    DM and symbolically for SX/MX."""
    Arr = c.ArrayInterface

    def arr_close(got, expected, msg):
      a = numpy.asarray(got, dtype=float)
      b = numpy.asarray(expected, dtype=float)
      self.assertEqual(a.shape, b.shape, "%s shape %s!=%s" % (msg, a.shape, b.shape))
      self.checkarray(ca.DM(a.reshape(-1, 1) if a.ndim else a),
                      ca.DM(b.reshape(-1, 1) if b.ndim else b), msg)

    x = Arr([[1, 2, 3], [4, 5, 6]])
    v = Arr([1, 2, 3])
    s = Arr(5)
    self.assertEqual((x.ndim, x.shape), (2, (2, 3)))
    self.assertEqual((v.ndim, v.shape), (1, (3,)))
    self.assertEqual((s.ndim, s.shape), (0, ()))

    # numpy-style indexing (rows, not columns)
    self.assertEqual(x[0].shape, (3,))
    arr_close(x[0], [1, 2, 3], "x[0] row")
    arr_close(x[1, 2], 6, "x[1,2] scalar")
    arr_close(x[:, 1], [2, 5], "x[:,1] column as 1-D")
    arr_close(v[1:3], [2, 3], "v slice")

    # iteration yields rows / elements
    rows = list(x)
    self.assertEqual(len(rows), 2)
    arr_close(rows[1], [4, 5, 6], "iter row1")
    self.assertEqual(len(list(v)), 3)

    # numpy right-aligned broadcasting (casadi alone refuses (1,3)+(2,3))
    arr_close(x + Arr([10, 20, 30]), [[11, 22, 33], [14, 25, 36]],
              "row-broadcast add")
    arr_close(x + 100, [[101, 102, 103], [104, 105, 106]], "scalar add")

    # matmul (numpy semantics across ndim)
    A = Arr([[1, 2], [3, 4]])
    arr_close(A @ Arr([1, 1]), [3, 7], "2D@1D")
    arr_close(Arr([1, 2, 3]) @ Arr([1, 1, 1]), 6, "1D@1D dot")

    # numpy dispatch returns wrapper instances
    r = numpy.sin(v)
    self.assertIsInstance(r, Arr)
    arr_close(r, numpy.sin([1, 2, 3]), "np.sin(1D)")
    arr_close(numpy.sum(x, axis=0), [5, 7, 9], "np.sum axis0")
    arr_close(numpy.transpose(x), [[1, 4], [2, 5], [3, 6]], "np.transpose")

    # setitem
    y = Arr([[0, 0, 0], [0, 0, 0]])
    y[0] = Arr([1, 2, 3])
    arr_close(y[0], [1, 2, 3], "row setitem")

    # symbolic: an ArrayInterface can be passed straight into a casadi API (the
    # typemap converts it -- a 1-D result maps to a column (n,1)); no
    # explicit to_casadi() needed.
    for cls in [ca.SX, ca.MX]:
      sym = cls.sym("M", 2, 3)
      ax = Arr(sym)
      expr = numpy.sin(ax[0])           # sin of the first row -> 1-D (3,)
      self.assertIsInstance(expr, Arr)
      f = ca.Function("f", [sym], [expr])  # <-- ArrayInterface straight in
      self.checkarray(f(ca.DM([[1, 2, 3], [4, 5, 6]])),
                      ca.DM(numpy.sin([1.0, 2.0, 3.0])),   # column (3,1)
                      "%s symbolic np.sin(row)" % cls.__name__)

  def test_experimental_ndarray(self):
    """issue #2959: dense N-D (>2-D) tensor support on the experimental
    array -- reshape, transpose(axes), broadcasting elementwise ops, axis
    reductions and stack/concatenate -- validated against numpy."""
    Arr = c.ArrayInterface
    np = numpy

    def same(label, got, ref):
      g = np.asarray(got, dtype=float)
      r = np.asarray(ref, dtype=float)
      self.assertEqual(g.shape, r.shape, "%s shape %s!=%s" % (label, g.shape, r.shape))
      self.assertTrue(np.allclose(g, r), "%s value" % label)

    A = np.arange(24.0).reshape(2, 3, 4)
    B = np.arange(24.0).reshape(2, 3, 4) + 100
    xA, xB = Arr(A), Arr(B)
    self.assertEqual((xA.ndim, xA.shape), (3, (2, 3, 4)))

    # reshape (incl. to/from >2-D and -1)
    same("reshape 3D->2D", np.reshape(xA, (6, 4)), np.reshape(A, (6, 4)))
    same("reshape 2D->3D", np.reshape(Arr(A.reshape(6, 4)), (2, 3, 4)), A)
    same("reshape -1", np.reshape(xA, (4, -1)), np.reshape(A, (4, -1)))

    # transpose with an axis permutation
    same("transpose default", np.transpose(xA), np.transpose(A))
    same("transpose (2,0,1)", np.transpose(xA, (2, 0, 1)), np.transpose(A, (2, 0, 1)))

    # elementwise + right-aligned broadcasting
    same("3D+3D", xA + xB, A + B)
    same("np.sin(3D)", np.sin(xA), np.sin(A))
    same("3D + (4,)", xA + Arr(np.arange(4.0)), A + np.arange(4.0))
    same("3D + (3,1)", xA + Arr(np.arange(3.0).reshape(3, 1)),
         A + np.arange(3.0).reshape(3, 1))

    # axis / multi-axis reductions
    same("sum axis0", np.sum(xA, axis=0), np.sum(A, axis=0))
    same("sum (0,2)", np.sum(xA, axis=(0, 2)), np.sum(A, axis=(0, 2)))
    same("mean (1,2)", np.mean(xA, axis=(1, 2)), np.mean(A, axis=(1, 2)))
    same("sum None", np.sum(xA), np.sum(A))

    # stack (new axis) / concatenate (existing axis)
    same("stack axis2", np.stack([xA, xB], axis=2), np.stack([A, B], axis=2))
    same("concat axis1", np.concatenate([xA, xB], axis=1), np.concatenate([A, B], axis=1))
    P = np.arange(6.0).reshape(2, 3)
    same("stack 2D->3D", np.stack([Arr(P), Arr(P + 9)]),
         np.stack([P, P + 9]))

    # symbolic dense N-D: reduce to <=2-D, then build a Function
    for cls in [ca.SX, ca.MX]:
      v = cls.sym("v", 24)
      t = np.reshape(Arr(v), (2, 3, 4))
      r = np.sum(t, axis=2)                 # (2, 3) -> back to casadi
      f = ca.Function("f", [v], [r.to_casadi()])  # pyright: ignore[reportAttributeAccessIssue]
      same("%s 3D sum axis2" % cls.__name__,
           f(ca.DM(A.reshape(-1, 1))), np.sum(A, axis=2))

    # an N-D array cannot be handed back to a 2-D casadi API
    w = ca.SX.sym("v", 24)
    with self.assertRaises(Exception):
      ca.Function("f", [w], [np.reshape(Arr(w), (2, 3, 4))])

  def test_issue4268(self):

    class Foo:
        def __rmatmul__(self, other):
            return 1
        def __radd__(self, other):
            return 2

    x= ca.MX.sym("x",2,2)

    f=Foo()
    print(x + f)
    print(x @ f)
    
    class Foo:
        def __req__(self, other):
            return 2
        def __eq__(self, other):  # pyright: ignore[reportIncompatibleMethodOverride]
            return 3
        def __rne__(self, other):
            return 4
        def __ne__(self, other):  # pyright: ignore[reportIncompatibleMethodOverride]
            return 5
              
    f = Foo()
    for X in [ca.SX,ca.MX]:
        x = X.sym("x")
        
        
        self.assertEqual(f==x,3)
        self.assertEqual(x==f,3)
        self.assertEqual(f!=x,5)
        self.assertEqual(x!=f,5)
    
    for x in [ca.SX.sym("x"),ca.SX.zeros(1,1),ca.SX.ones(1,1)]:
        for y in [ca.MX.sym("y"),ca.MX.zeros(1,1),ca.MX.ones(1,1)]:
        
            with self.assertInException("Cannot compare SX and MX objects for equality"):
                x==y  # pyright: ignore[reportOperatorIssue]
            with self.assertInException("Cannot compare SX and MX objects for equality"):
                y==x  # pyright: ignore[reportOperatorIssue]
            with self.assertInException("Cannot compare SX and MX objects for inequality"):
                x!=y  # pyright: ignore[reportOperatorIssue]
            with self.assertInException("Cannot compare SX and MX objects for inequality"):
                y!=x  # pyright: ignore[reportOperatorIssue]

  def test_conversion_operators(self):
    self.message("COnversion operations")


    def doit(z,s,fun):
      if type(z) in [type(ca.SX()),type(ca.SX())]:
        ztype = [type(ca.SX()),type(ca.SX())]

      if type(z) in [type(ca.MX())]:
        ztype = [type(ca.MX())]

      r = fun(z,s)

      if type(z) is type(ca.SX()) and type(s) is type(ca.SX()):
        self.assertTrue(type(r) is type(ca.SX()))


      self.assertTrue(type(r) in ztype,"Expected %s but got %s" % (str(ztype),str(type(r))))

    def tests(z,s):
      doit(z,s,lambda z,s: s>=z)
      doit(z,s,lambda z,s: s>z)
      doit(z,s,lambda z,s: s<=z)
      doit(z,s,lambda z,s: s<z)
      doit(z,s,lambda z,s: s==z)
      doit(z,s,lambda z,s: s!=z)

    nums = [np.array([[1,2],[3,4]]),ca.DM([[1,2],[3,4]]), ca.DM(4), np.array(4),4.0,4]

    ## numeric & SX
    for s in nums:
      for z in [ca.SX.sym("x"), ca.SX.sym("x"), ca.SX.sym("x",2,2)]:
        print("z = %s, s = %s" % (str(z),str(s)))
        print("  z = %s, s = %s" % (type(z),type(s)))
        tests(z,s)

    # numeric & MX
    for s in nums:
      for z in [ca.MX.sym("x",2,2)]:
        print("z = %s, s = %s" % (str(z),str(s)))
        print("  z = %s, s = %s" % (type(z),type(s)))
        tests(z,s)

    # SX & SX
    for s in [ca.SX.sym("x"), ca.SX.sym("x"), ca.SX.sym("x",2,2)]:
      for z in [ca.SX.sym("x"),ca.SX.sym("x"), ca.SX.sym("x",2,2)]:
        print("z = %s, s = %s" % (str(z),str(s)))
        print("  z = %s, s = %s" % (type(z),type(s)))
        tests(z,s)

    # MX & MX
    for s in [ca.MX.sym("x"),ca.MX.sym("x",2,2)]:
      for z in [ca.MX.sym("x"),ca.MX.sym("x",2,2)]:
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
    w=ca.DM(goal)
    self.checkarray(w,goal,"Constructor")

    for name, value in list(test.items()):
      w.nz[:] = value
      self.checkarray(w,goal,"name")

    test={
      "array2ddouble" : np.array([goallist],dtype=double).T,
      "array2dint" : np.array([goallist]).T,
    }
    w=ca.DM(goal)
    self.checkarray(w,goal,"Constructor")

    for name, value in list(test.items()):
      w[:,:] = value
      self.checkarray(w,goal,"name")

  def test_DMSXcast(self):
    self.message("Casting DM to SX")

    W = ca.SX(ca.DM([[1,2,3],[4,5,6]]))

    self.assertEqual(W.size1(),2)
    self.assertEqual(W.size2(),3)

  def test_DMMXcast(self):
    self.message("Casting DM to MX")
    W = ca.MX(ca.DM([[1,2,3],[4,5,6]]))

    self.assertEqual(W.size1(),2)
    self.assertEqual(W.size2(),3)

  def test_DMSX(self):
    self.message("Casting DM to SX")

    w = ca.DM([[1,2,3],[4,5,6]])
    x = ca.SX.sym("x")

    f = ca.Function("f", [x],[w])

    W = f(*f.sx_in())
    self.assertEqual(W.size1(),2)
    self.assertEqual(W.size2(),3)

  def test_DMMX(self):
    self.message("Casting DM to MX")
    w = ca.DM([[1,2,3],[4,5,6]])
    x = ca.MX.sym("x")

    f = ca.Function("f", [x],[w])

    W = f(*f.mx_in())

    self.assertEqual(W.size1(),2)
    self.assertEqual(W.size2(),3)

  def test_setgetslice(self):
    self.message("set/get on DM using slices")

    w = ca.DM([[0,0]])

    A = np.array([[1.0,2],[3,4]])
    B = np.array([[4.0,5],[6,7]])

    w[:,:] = A[0,:]
    self.checkarray(w,A[0,:].reshape((1,-1)),"set")
    B[:1,:] = w
    self.checkarray(B[0,:],A[0,:],"get")

    w = ca.DM([[0],[0]])


    w[:,:] = A[:,0]
    self.checkarray(w,A[:,0],"set")
    B[:,:1] = w
    self.checkarray(B[:,0],A[:,0],"get")

    w = ca.DM([[1,2],[3,4]])
    A = np.zeros((8,7))
    B = np.zeros((8,7))
    A[2:7:3,:7:4] = w
    B[2:7:3,:7:4] = w
    self.checkarray(A,B,"get")

  def test_vertcatprecedence(self):
    self.message("Argument precedence DM")
    a = ca.DM([1,2])
    self.assertTrue(isinstance(ca.vertcat(*[a,a]),ca.DM))

    a = ca.DM([1,2])
    self.assertTrue(isinstance(ca.vertcat(*[a,[1,2,3]]),ca.DM))


    a = ca.MX([1,2])
    print(ca.vertcat(*[a,[1,2,3]]))
    self.assertTrue(isinstance(ca.vertcat(*[a,[1,2,3]]),ca.MX))

  def test_issue190(self):
    self.message("regression test issue #190")
    x=ca.SX.sym("x")
    x * np.array(1)
    x * np.array(1.2)

    ca.SX.sym("x") * np.array(1.0)
    ca.MX.sym("x") * np.array(1.0)

  def test_array_cat(self):
    ca.horzcat(*(ca.SX.sym("x",4,3),np.ones((4,3))))


  def test_issue(self):
    self.message("std::vector<double> typemap.")
    a = np.array([0,2,2,3])
    b = np.array([0.738,0.39,0.99])
    ca.DM(ca.Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.39,0.99])
    ca.DM(ca.Sparsity(4,3,(0,2,2,3),[1,2,1]),[0.738,0.39,0.99])
    ca.DM(ca.Sparsity(4,3,list(a),[1,2,1]),[0.738,0.39,0.99])
    ca.DM(ca.Sparsity(4,3,a,[1,2,1]),[0.738,0.39,0.99])
    ca.DM(ca.Sparsity(4,3,[0,2,2,3],[1,2,1]),(0.738,0.39,0.99))
    ca.DM(ca.Sparsity(4,3,[0,2,2,3],[1,2,1]),list(b))
    ca.DM(ca.Sparsity(4,3,[0,2,2,3],[1,2,1]),b)

  def test_issue314(self):
    self.message("regression test for #314: SX sparsity constructor")
    ca.SX(ca.Sparsity.diag(3),[1,2,3])
  def test_setAll_365(self):
    self.message("ticket #365: DMAtrix.setAll does not work for 1x1 Matrices as input")
    for i in [ca.DM(4),np.array([4]),np.array(4),np.array([[4]]),4,4.0]:
      m = ca.DM.ones(5,5)
      m[:,:] = i
      self.checkarray(m,ca.DM.ones(5,5)*4)

  @unittest.skipIf(sys.version_info>=(3,0), "To lazy to fix")
  def test_issue570(self):
    self.message("Issue #570: long int")
    longint = 10**50
    print(type(longint))
    print(ca.SX.sym('x') + longint)
    print(longint + ca.SX.sym('x'))
    print(ca.SX.sym('x') + longint)
    print(longint + ca.SX.sym('x'))

  def test_casting_DM(self):
    self.message("casting DM")

    print("cast A")
    x = ca.SX.sym("x")
    f = ca.Function("f", [x],[x])
    class Foo:
      def __DM__(self):
        return ca.DM([4])

    self.assertEqual(f(Foo()),4)
    print("cast B")

    class Foo:
      def __DM__(self):
        return ca.SX([4])


    self.assertRaises(TypeError if systemswig else NotImplementedError,lambda :f(Foo()))  # pyright: ignore[reportCallIssue,reportArgumentType]
    print("cast C")

    class Foo:
      def __DM__(self):
        raise Exception("15")

    self.assertRaises(TypeError if systemswig else NotImplementedError,lambda :f(Foo()))  # pyright: ignore[reportCallIssue,reportArgumentType]
    print("cast D")

    class Foo:
      pass

    self.assertRaises(TypeError if systemswig else NotImplementedError,lambda :f(Foo()))  # pyright: ignore[reportCallIssue,reportArgumentType]
    print("cast E")

  def test_casting_SX(self):
    self.message("casting SX")


    x = ca.SX.sym("x")

    class Foo:
      def __SX__(self):
        return x

    ca.Function("tmp", [x],[Foo()])

    class Foo:
      def __SX__(self):
        return ca.MX.sym("x")

    self.assertRaises(TypeError if systemswig else NotImplementedError,lambda : ca.Function("tmp", [x],[Foo()]))  # pyright: ignore[reportCallIssue,reportArgumentType]

    class Foo:
      def __SX__(self):
        raise Exception("15")

    self.assertRaises(TypeError if systemswig else NotImplementedError,lambda : ca.Function("tmp", [x],[Foo()]))  # pyright: ignore[reportCallIssue,reportArgumentType]

    class Foo:
      pass

    self.assertRaises(TypeError if systemswig else NotImplementedError,lambda :ca.Function("tmp", [x],[Foo()]))  # pyright: ignore[reportCallIssue,reportArgumentType]


  def test_casting_MX(self):
    self.message("casting MX")


    x = ca.MX.sym("x")

    class Foo:
      def __MX__(self):
        return x

    ca.Function("tmp", [x],[Foo()])

    class Foo:
      def __MX__(self):
        return ca.SX.sym("x")

    with self.assertRaises(TypeError if systemswig else NotImplementedError):
      ca.Function("tmp", [x],[Foo()])  # pyright: ignore[reportCallIssue,reportArgumentType]

    class Foo:
      def __MX__(self):
        raise Exception("15")

    with self.assertRaises(Exception):
      ca.Function("tmp", [x],[Foo()])  # pyright: ignore[reportCallIssue,reportArgumentType]

    class Foo:
      pass

    with self.assertRaises(TypeError if systemswig else NotImplementedError):
      ca.Function("tmp", [x],[Foo()])  # pyright: ignore[reportCallIssue,reportArgumentType]

  def test_OUTPUT(self):
    self.message("OUTPUT typemap")
    a = ca.SX.sym("A",3,3)
    
    self.assertTrue(isinstance(ca.qr(a),tuple([tuple]+([list] if swig4 else []))))  # pyright: ignore[reportArgumentType]

  def test_cvar(self):
    self.message("We must not have cvar, to avoid bug #652")
    # Wrap all static global things in #ifdef SWIG
    with self.assertRaises(Exception):
      cvar  # pyright: ignore[reportUndefinedVariable]

  def test_ufuncsum(self):
    self.message("ufunc.add")

    self.checkarray(ca.DM(np.sum(ca.DM([1,2,3]).nonzeros())),ca.DM(6))

  def test_sxmatrix(self):

    def val(a):
      f = ca.Function("f", [],[a])
      f_out = f.call([])
      return f_out[0]

    for i in [ca.SX(1),1,1.0]:
      a = np.array([[1,2],[3,4]])
      print(val(ca.SX(a)))
      print(val(ca.SX(a.T)))

      self.checkarray(val(ca.SX(a)),ca.DM([[1,2],[3,4]]))
      self.checkarray(val(ca.SX(a.T).T),ca.DM([[1,2],[3,4]]))

      if check_matrix:
        a = numpy.matrix([[1,2],[3,4]])

        print(val(ca.SX(a)))
        print(ca.DM([[1,2],[3,4]]))

        self.checkarray(val(ca.SX(a)),ca.DM([[1,2],[3,4]]))
        self.checkarray(val(ca.SX(a.T).T),ca.DM([[1,2],[3,4]]))

  def test_issue1158(self):
    A = numpy.zeros((0,2))
    a = ca.DM(A)
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
          csc_matrix(([1.0,3.0,2.0,4.0],[0,1,0,1],[0,2,4]),shape=(2,2),dtype=numpy.double),  # pyright: ignore[reportUnboundVariable]
          csc_matrix(([1,3,2,4],[0,1,0,1],[0,2,4]),shape=(2,2),dtype=numpy.intc),  # pyright: ignore[reportUnboundVariable]
          ca.DM([[1,2],[3,4]]).sparse()
      ]


    for D in Ds:
      print(D)
      d = ca.DM.ones(2,2)

      x = ca.SX.sym("x",d.sparsity())
      f = ca.Function("f", [x],[x])
      fin = ca.DM(d.sparsity(),0)
      fin[:,:] = D

      self.checkarray(fin,ca.DM([[1,2],[3,4]]))
      d[:,:] = D
      self.checkarray(d,ca.DM([[1,2],[3,4]]))

  def test_issue1217(self):
    a = np.array([0,ca.SX.sym("x")])

    print(ca.if_else(0,a,a))

  def test_issue1373(self):
    print(np.array(ca.DM([2])))
    print(np.array(ca.DM([1,2,3.0])))

  def test_None(self):
    #self.assertFalse(None==DM(3))
    b = np.atleast_2d(None)  # pyright: ignore[reportCallIssue,reportArgumentType]
    with self.assertRaises(TypeError if systemswig else NotImplementedError):
      c = ca.repmat(b, 1, 1)

  @requires_nlpsol("ipopt")
  def testGenericTypeBoolean(self):
    x=ca.SX.sym("x")
    with self.assertRaises(RuntimeError):
      solver = ca.nlpsol("mysolver", "ipopt", {"x":x,"f":x**2}, {"ipopt": {"acceptable_tol": ca.SX.sym("x")}})

    ca.nlpsol("mysolver", "ipopt", {"x":x,"f":x**2}, {"ipopt": {"acceptable_tol": 1}})
  
  """
  def test_to_longlong(self):
    a = IM(10)


    b = a**15

    self.assertEqual(int(b),10**15)
  """
  
  def test_issue_3214(self):
    import numpy as np
    import casadi as ca
    # `numpy_array - DM` is operator arithmetic (subtract ufunc): it stays
    # casadi-typed with casadi's strict shape rules, regardless of
    # setNumpyMode.  A (2,) minus a (1,2) is a dimension mismatch.
    with self.assertRaises(RuntimeError):
        print(np.array([-1.0, 0.0]) - ca.DM([[0.1, 0.0]]))

    r = np.array([-1.0, 0.0]) - ca.DM([0.1, 0.0])
    self.checkarray(r,ca.DM([-1.1,0]))

  def test_buglonglong(self):
    x = ca.SX.sym("x")

    ca.jacobian(x/1.458151064450277e-12,x)

  def test_issue273(self):
    state = ca.SX.sym('H_body_world', 6, 6)

    H_body_world = np.array([[ca.cos(state[2, 0])*ca.cos(state[1, 0]), -ca.sin(state[2, 0])*ca.cos(state[1, 0]), ca.sin(state[1, 0]), state[3, 0]],
                             [ca.sin(state[0, 0])*ca.sin(state[1, 0])*ca.cos(state[2, 0]) + ca.sin(state[2, 0])*ca.cos(state[0, 0]), -ca.sin(state[0, 0])*ca.sin(state[2, 0])*ca.sin(state[1, 0]) + ca.cos(state[0, 0])*ca.cos(state[2, 0]), -ca.sin(state[0, 0])*ca.cos(state[1, 0]), state[4, 0]], 
                             [ca.sin(state[0, 0])*ca.sin(state[2, 0]) - ca.sin(state[1, 0])*ca.cos(state[0, 0])*ca.cos(state[2, 0]), ca.sin(state[0, 0])*ca.cos(state[2, 0]) + ca.sin(state[2, 0])*ca.sin(state[1, 0])*ca.cos(state[0, 0]), ca.cos(state[0, 0])*ca.cos(state[1, 0]), state[5, 0]],
                             [0, 0, 0, 1]])

    print(H_body_world)
    print(ca.SX(H_body_world))

  def test_issue_2625(self):
    # The casadi-typed inner/outer/reduction semantics below are the
    # type-preserving bridge (issue #2959); opt in so numeric DM inputs
    # take the casadi handler instead of densifying to numpy.
    ca.GlobalOptions.setNumpyMode(1)
    try:
      self._body_issue_2625()
    finally:
      ca.GlobalOptions.setNumpyMode(0)

  def _body_issue_2625(self):
    # Under setNumpyMode the numpy ops return an
    # casadi.ArrayInterface; checkarray densifies it via __array__, while
    # casadi APIs need the underlying value (to_casadi()).
    Arr = c.ArrayInterface
    def _uw(e):
      return e.to_casadi() if isinstance(e, Arr) else e

    # np.inner is a true inner product, not the outer-product-shaped result the
    # old buggy fallback produced.  np.outer is the outer product as expected.
    self.checkarray(np.inner(ca.DM([1,0,1]),ca.DM([1,0,1])), np.array(2))
    self.checkarray(np.outer(ca.DM([1,0,1]),ca.DM([1,0,1])), np.array([[1,0,1],[0,0,0],[1,0,1]]))

    # Reductions on DM produce shape-correct (axis-dropped) numpy values.
    self.checkarray(np.cumsum(ca.DM([1,0,1])), np.array([1,1,2]))
    self.checkarray(np.sum(ca.DM([1,0,1])),np.array(2))
    self.assertFalse(np.all(ca.DM([1,0,1])))
    self.assertTrue(np.any(ca.DM([1,0,1])))
    self.checkarray(np.sum(ca.DM([[1,0,1],[0,1,0]])),np.array(3))
    self.assertFalse(np.all(ca.DM([[1,0,1],[0,1,0]])))
    self.assertTrue(np.any(ca.DM([[1,0,1],[0,1,0]])))

    # Reductions and inner/outer now work natively on symbolic types too.
    def _eval(expr):
      return ca.Function("f", [], [_uw(expr)])()["o0"]

    for M in [ca.SX, ca.MX]:
      v = M([1, 0, 1])
      self.checkarray(_eval(np.inner(v, v)), ca.DM(2), digits=14)
      self.checkarray(_eval(np.outer(v, v)),
                      ca.DM([[1, 0, 1], [0, 0, 0], [1, 0, 1]]), digits=14)
      self.checkarray(_eval(ca.vec(_uw(np.cumsum(v)))), ca.DM([1, 1, 2]), digits=14)
      self.checkarray(_eval(np.sum(v)),                     ca.DM(2),         digits=14)

    # np.all / np.any are built from mmin/mmax of (x != 0), so they work on
    # MX too (no logic_all/logic_any core primitive needed).
    for M in [ca.SX, ca.MX]:
      self.checkarray(_eval(np.all(M([1, 1, 1]))), ca.DM(1))
      self.checkarray(_eval(np.any(M([0, 0, 1]))), ca.DM(1))
      self.checkarray(_eval(np.all(M([1, 0, 1]))), ca.DM(0))

if __name__ == '__main__':
    unittest.main()
