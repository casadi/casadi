# coding=utf-8
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
import unittest
from types import *
from helpers import *
import pickle
from operator import itemgetter
import sys
from casadi.tools import capture_stdout

scipy_available = True
try:
    import scipy.special
    from scipy.linalg import expm
except:
    scipy_available = False

class Misctests(casadiTestCase):

  def test_issue179B(self):
    self.message('Regression test #179 (B)')
    def calc_sparsity():
      x = ca.SX.sym("x")
      f = ca.Function('f', [x], [x ** 2])
      return f.jac_sparsity(0, 0)

    def print_sparsity():
        sparsity = calc_sparsity()
        str(sparsity) # Segfault

    print_sparsity()

  def test_issue4216(self):
    self.message('Regression test #4216: complex numpy inputs should raise, not segfault')
    self.assertRaises(TypeError, lambda: 0.5j * ca.MX(1))  # pyright: ignore[reportCallIssue,reportArgumentType,reportOperatorIssue,reportAttributeAccessIssue]

    c = numpy.array([3.+2j])
    self.assertRaises((TypeError, NotImplementedError), lambda: ca.DM(c))  # pyright: ignore[reportCallIssue,reportArgumentType,reportOperatorIssue,reportAttributeAccessIssue]
    self.assertRaises((TypeError, NotImplementedError), lambda: ca.SX(c))  # pyright: ignore[reportCallIssue,reportArgumentType,reportOperatorIssue,reportAttributeAccessIssue]
    self.assertRaises((TypeError, NotImplementedError), lambda: ca.MX(c))  # pyright: ignore[reportCallIssue,reportArgumentType,reportOperatorIssue,reportAttributeAccessIssue]

    cv = 0.5j * numpy.ones(2)
    # MX has no implicit numeric conversion so matmul reaches the casadi
    # path and must raise TypeError. For DM/SX, numpy's __array__ path
    # may instead raise ValueError on shape mismatch; either is fine here
    # as long as no segfault occurs.
    self.assertRaises((TypeError, NotImplementedError), lambda: cv @  ca.MX(2))  # pyright: ignore[reportCallIssue,reportArgumentType,reportOperatorIssue,reportAttributeAccessIssue]
    for d in (ca.DM(2), ca.SX(2)):
      self.assertRaises((TypeError, NotImplementedError, ValueError),
                        lambda d=d: cv @ d)  # pyright: ignore[reportCallIssue,reportArgumentType,reportOperatorIssue,reportAttributeAccessIssue]

  def test_sanity(self):
    ca.DM(ca.Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.39,0.99])
    self.assertRaises(RuntimeError,lambda : ca.DM(ca.Sparsity(4,4,[0,2,2,3],[1,2,1]),[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : ca.DM(ca.Sparsity(4,3,[0,2,2,12],[1,2,1]),[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : ca.DM(ca.Sparsity(4,3,[-10,2,2,3],[1,2,1]),[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : ca.DM(ca.Sparsity(4,3,[0,2,2,3],[8,2,1]),[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : ca.DM(ca.Sparsity(4,3,[0,2,2,3],[-3,2,1]),[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : ca.DM(ca.Sparsity(4,3,[0,2,2,3],[1,2,1,2]),[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : ca.DM(ca.Sparsity(4,3,[0,2,0,3],[1,2,1]),[0.738,0.39,0.99]))

  def test_copyconstr_norefcount(self):
    self.message("Copy constructor for non-refcounted classes")
    x = ca.DM.ones(2,3)

    y = ca.DM(x)
    x[0,0] = 5

    self.assertFalse(id(x)==id(y))
    self.assertEqual(x[0,0],5)
    self.assertEqual(y[0,0],1)

  def test_copyconstr_refcount(self):
    self.message("Copy constructor for refcounted classes")
    x = ca.Sparsity.diag(4)

    y = ca.Sparsity(x)

    x.resize(2,8)

    self.assertFalse(id(x)==id(y))
    self.assertTrue(x.numel(),y.numel())
    self.checkarray(x.shape,(2,8),"shape")
    self.checkarray(y.shape,(4,4),"shape")

  def test_copy_norefcount(self):
    self.message("Shallow copy for non-refcounted classes")
    import copy

    x = ca.DM.ones(2,3)

    y = copy.copy(x)
    x[0,0] = 5

    self.assertFalse(id(x)==id(y))
    self.assertEqual(x[0,0],5)
    self.assertEqual(y[0,0],1)

  def test_copy_refcount(self):
    self.message("Shallow copy for refcounted classes")
    import copy
    x = ca.Sparsity.diag(4)

    y = copy.copy(x)

    x.resize(2,8)

    self.assertFalse(id(x)==id(y))
    self.assertTrue(x.numel(),y.numel())
    self.checkarray(x.shape,(2,8),"shape")
    self.checkarray(y.shape,(4,4),"shape")

  def test_deepcopy_norefcount(self):
    self.message("Deep copy for non-refcounted classes")
    import copy

    x = ca.DM.ones(2,3)

    y = copy.deepcopy(x)
    x[0,0] = 5

    self.assertFalse(id(x)==id(y))
    self.assertEqual(x[0,0],5)
    self.assertEqual(y[0,0],1)

  def test_deepcopy_refcount(self):
    self.message("Deep copy for refcounted classes")
    import copy
    x = ca.Sparsity.diag(4)

    y = copy.deepcopy(x)

    x.resize(2,8)

    self.assertFalse(id(x)==id(y))
    self.assertTrue(x.numel(),y.numel())
    self.checkarray(x.shape,(2,8),"shape")
    self.checkarray(y.shape,(4,4),"shape")

  @requiresPlugin(ca.nlpsol,"ipopt")
  def test_options_introspection(self):
    self.message("options introspection")
    x=ca.SX.sym("x")
    nlp = {'x':x, 'f':x**2}
    i = ca.nlpsol('i', "ipopt", nlp)

    opts = i.optionNames()  # pyright: ignore[reportAttributeAccessIssue]
    self.assertTrue(isinstance(opts,list))

    n = opts[0]
    self.assertTrue(type(n)==type(""))

    n = "monitor"

    d = i.optionDescription(n)  # pyright: ignore[reportAttributeAccessIssue]
    self.assertTrue(type(d)==type(""))
    self.assertTrue(not("d"=="N/A"))

    d = i.optionTypeName(n)  # pyright: ignore[reportAttributeAccessIssue]
    self.assertEqual(d,"OT_STRINGVECTOR")

    #d = i.optionAllowed(n)

  def test_pickling(self):

    a = ca.Sparsity.lower(4)
    s = pickle.dumps(a)
    b = pickle.loads(s)
    self.assertTrue(a==b)

    a = ca.Sparsity()
    s = pickle.dumps(a)
    b = pickle.loads(s)
    self.assertTrue(a.is_null())

    a = ca.DM(ca.Sparsity.lower(4),list(range(10)))
    s = pickle.dumps(a)
    b = pickle.loads(s)
    self.checkarray(a,b)

  def test_exceptions(self):
    if systemswig: return

    try:
      ca.nlpsol(123)  # pyright: ignore[reportCallIssue,reportArgumentType,reportOperatorIssue,reportAttributeAccessIssue]
      self.assertTrue(False)
    except NotImplementedError as e:
      e_message = str(e)
      assert "nlpsol(str,str,dict:MX,dict)" in e_message
      assert "You have: '(int)'" in e_message
      assert "::" not in e_message
      assert "std" not in e_message

    try:
      ca.vcat(123)  # pyright: ignore[reportCallIssue,reportArgumentType,reportOperatorIssue,reportAttributeAccessIssue]
      self.assertTrue(False)
    except NotImplementedError as e:
      e_message = str(e)
      assert "vertcat([SX]" in e_message
      assert "vertcat([DM]" in e_message
      assert "You have: '(int)'" in e_message
      assert "::" not in e_message
      assert "std" not in e_message

    try:
      ca.substitute(123)  # pyright: ignore[reportCallIssue,reportArgumentType,reportOperatorIssue,reportAttributeAccessIssue]
      self.assertTrue(False)
    except NotImplementedError as e:
      e_message = str(e)
      assert "substitute(SX,SX,SX)" in e_message
      assert "substitute([SX],[SX],[SX])" in e_message
      assert "You have: '(int)'" in e_message
      assert "::" not in e_message
      assert "std" not in e_message

    try:
      ca.load_nlpsol(132)  # pyright: ignore[reportCallIssue,reportArgumentType,reportOperatorIssue,reportAttributeAccessIssue]
      self.assertTrue(False)
    except NotImplementedError as e:
      e_message = str(e)
      assert "load_nlpsol(str)" in e_message
      assert "::" not in e_message
      assert "std" not in e_message

    x=ca.SX.sym("x")

    try:
      [x]+ x  # pyright: ignore[reportCallIssue,reportArgumentType,reportOperatorIssue,reportAttributeAccessIssue]
      self.assertTrue(False)
    except TypeError as e:
      e_message = str(e)

    try:
      x + [x]  # pyright: ignore[reportCallIssue,reportArgumentType,reportOperatorIssue,reportAttributeAccessIssue]
      self.assertTrue(False)
    except TypeError as e:
      e_message = str(e)

    try:
      x.reshape(2)
      self.assertTrue(False)
    except NotImplementedError as e:
      e_message = str(e)
      assert "reshape(SX,(int,int))" in e_message

    try:
      x.reshape(("a",2))  # pyright: ignore[reportCallIssue,reportArgumentType,reportOperatorIssue,reportAttributeAccessIssue]
      self.assertTrue(False)
    except NotImplementedError as e:
      e_message = str(e)
      assert "You have: '(SX,(str,int))'" in e_message

    try:
      ca.diagsplit("s")  # pyright: ignore[reportCallIssue,reportArgumentType,reportOperatorIssue,reportAttributeAccessIssue]
      self.assertTrue(False)
    except NotImplementedError as e:
      e_message = str(e)
      assert "diagsplit(SX,int)" in e_message
      assert "diagsplit(DM,int)" in e_message

    try:
      ca.DM("df")  # pyright: ignore[reportCallIssue,reportArgumentType,reportOperatorIssue,reportAttributeAccessIssue]
      self.assertTrue(False)
    except NotImplementedError as e:
      e_message = str(e)
      assert "  DM(" in e_message

    try:
      ca.vertcat(1,ca.SX.sym('x'),ca.MX.sym('x'))  # pyright: ignore[reportCallIssue,reportArgumentType,reportOperatorIssue,reportAttributeAccessIssue]
      self.assertTrue(False)
    except NotImplementedError as e:
      e_message = str(e)
      assert "vertcat(" in e_message

  def test_getscheme(self):
    x = ca.SX.sym("x")
    p = ca.SX.sym("p")

    F = ca.Function('F', [x, p], [x+p, x**2], ['x', 'p'], ['f', 'g'])

    fc = F(x=3,p=4)
    f = fc['f']
    self.checkarray(f,ca.DM([7]))
    g = fc['g']
    self.checkarray(g,ca.DM([9]))
    [f,g] = itemgetter('f','g')(fc)
    self.checkarray(f,ca.DM([7]))
    self.checkarray(g,ca.DM([9]))
    [g,f] = itemgetter('g','f')(fc)
    self.checkarray(f,ca.DM([7]))
    self.checkarray(g,ca.DM([9]))

  def test_assertions(self):

    x = ca.MX.sym("x")

    z = x**2

    z = z.attachAssert(z>3,"x must be larger than 3")

    v = ca.sin(z)

    f = ca.Function('f', [x],[v])

    print(f)

    f_out = f(-6)
    try :
      f_out = f(1)
    except Exception as e:
      print(str(e))
      self.assertTrue("x must be larger than 3" in str(e))

  def test_doc(self):
    doc = ca.nlpsol.__doc__ or ""
    self.assertTrue("jac_penalty" in doc) # FunctionInternal
    self.assertTrue("print_time" in doc)  # ProtoFunction
    self.assertTrue( doc.count("print_time")==1)  # ProtoFunction

  @requires_nlpsol("ipopt")
  def test_output(self):
    with capture_stdout() as result:
      ca.DM([1,2]).print_dense()


    assert "2" in result[0]

    x=ca.SX.sym("x")
    f = {'x':x, 'f':x**2}
    solver = ca.nlpsol("solver", "ipopt",f)
    with capture_stdout() as result:
      solver_out = solver(x0=0)

    assert "Number of nonzeros in equality constraint" in result[0]
    assert "iter    objective    inf_pr" in result[0]

    with capture_stdout() as result:
      try:
        solver = ca.nlpsol("solver","foo",f)
      except:
        pass

    assert "casadi_nlpsol_foo" in result[1]

  def test_serialize(self):

    x = ca.Sparsity.upper(3)

    y = ca.SX.sym("y") # nested

    si = ca.StringSerializer()
    si.pack([x])

    z = y
    for i in range(10000):
      z = ca.sin(z)
    e = ca.vertcat(y,z,2*z)
    fref = ca.Function('f',[y],[e])
    si.pack(e)
    data = si.encode()
    si.pack([x])

    si = ca.StringDeserializer(data)

    spx = si.unpack()
    self.assertTrue(spx[0]==x)
    A = si.unpack()
    f = ca.Function('f',[A[0]],[A])
    self.checkfunction_light(f,fref,[7])
    with self.assertInException("end of stream"):
      si.unpack()
    with self.assertInException("end of stream"):
      si.unpack()

    si = ca.FileSerializer("foo.dat")
    si.pack([x])
    z = ca.sin(y)
    e = ca.vertcat(y,z,2*z)
    fref = ca.Function('f',[y],[e])
    si.pack(e)
    si = None
    si = ca.FileDeserializer("foo.dat")

    spx = si.unpack()
    self.assertTrue(spx[0]==x)
    A = si.unpack()
    f = ca.Function('f',[A[0]],[A])
    self.checkfunction_light(f,fref,[7])
    with self.assertInException("end of stream"):
      si.unpack()
    with self.assertInException("end of stream"):
      si.unpack()


    x = ca.SX.sym('x')
    s = ca.StringSerializer()
    s.pack(x)
    data1 = s.encode()
    s.pack(ca.sin(x))
    data2 = s.encode()

    s = ca.StringDeserializer(data1)
    a = s.unpack()
    with self.assertInException("end of stream"):
      s.unpack()
    s.decode(data2)
    b = s.unpack()
    with self.assertInException("end of stream"):
      s.unpack()
    si = None

    si = ca.FileSerializer("foo.dat")
    f = ca.Function("f",[x],[x**2])
    g = ca.Function("g",[x],[x**3])
    si.pack([f,g])
    si = None
    with self.assertInException("File is not loadable with 'load'. Use 'FileDeserializer' instead."):
      r = ca.Function.load("foo.dat")

    si = ca.FileDeserializer("foo.dat")
    print(si.unpack())

    si = None

    si = ca.FileSerializer("foo.dat")
    f = ca.Function("f",[x],[x**2])
    si.pack(f)
    si = None
    r = ca.Function.load("foo.dat")
    print(r)

    f.save("foo.dat")
    si = ca.FileDeserializer("foo.dat")
    print(si.unpack())

  def test_print_time(self):


    x = ca.MX.sym("x")
    f = x**2


    for print_time in [True,False]:
      included = ["t_wall"] if print_time else []
      excluded = [] if print_time else ["t_wall"]
      opts = {"print_time":print_time}
      ff = ca.Function("f",[x],[f],opts)
      with self.assertOutput(included, excluded):
        ff(3)

      if ca.has_nlpsol("ipopt"):
        solver = ca.nlpsol("nlpsol","ipopt",{"x":x,"f":f},opts)
        with self.assertOutput(included, excluded):
          solver()


      if ca.has_conic("ipopt"):
        solver = ca.qpsol("qpsol","qpoases",{"x":x,"f":f},opts)
        with self.assertOutput(included, excluded):
          solver()

      solver = ca.rootfinder("rootfinder","newton",{"x":x,"g":x},opts)
      with self.assertOutput(included, excluded):
        solver()

      solver = ca.integrator("integrator","rk",{"x":x,"ode":f},opts)
      with self.assertOutput(included, excluded):
        solver()


      integr_options = dict(opts)
      integr_options["simplify"] = True
      solver = ca.integrator("integrator","rk",{"x":x,"ode":f},integr_options)
      with self.assertOutput(included, excluded):
        solver()

      A = ca.DM.rand(3,3)
      if ca.has_linsol("lapacklu"):
        with self.assertOutput(included, excluded):
          ca.solve(A,ca.vertcat(1,2,3),"lapacklu",opts)


      Amx = ca.MX(ca.DM.rand(3,3))
      if ca.has_linsol("lapacklu"):
        with self.assertOutput(included, excluded):
          ca.evalf(ca.solve(Amx,ca.vertcat(1,2,3),"lapacklu",opts))

      if ca.has_linsol("lapacklu"):
        solver = ca.Linsol("linsol","lapacklu",A.sparsity(),opts)
        with self.assertOutput(included, excluded):
          solver.solve(A,ca.vertcat(1,2,3))

  def test_error_formatting(self):
    x = ca.MX.sym("x")

    f = ca.Function('f',[x],[x.attachAssert(0,"Hey \\W %d")])

    with self.assertInException("Hey \\W %d"):
        f()

  @memory_heavy()
  def test_record_time(self):


    x = ca.MX.sym("x")
    f = x**2
    for i in range(10000):
      f =ca.sin(f)*f
    f=f+x**2

    opts = {"record_time":True}
    ff = ca.Function("f",[x],[f],opts)
    ff(3)
    self.assertTrue("t_proc_total" in ff.stats())
    self.assertTrue(ff.stats()["t_wall_total"]>=0)

    if ca.has_nlpsol("ipopt"):
      solver = ca.nlpsol("nlpsol","ipopt",{"x":x,"f":f},opts)
      solver()
      self.assertTrue("t_proc_total" in solver.stats())
      self.assertTrue(solver.stats()["t_wall_total"]>0)

    if ca.has_conic("qpoases"):
      solver = ca.qpsol("qpsol","qpoases",{"x":x,"f":f},opts)
      solver()
      self.assertTrue("t_proc_total" in solver.stats())
      self.assertTrue(solver.stats()["t_wall_total"]>0)

    solver = ca.rootfinder("rootfinder","newton",{"x":x,"g":x},opts)
    solver()
    self.assertTrue("t_proc_total" in solver.stats())
    self.assertTrue(solver.stats()["t_wall_total"]>=0)

    solver = ca.integrator("integrator","rk",{"x":x,"ode":f},opts)
    solver()
    self.assertTrue("t_proc_total" in solver.stats())
    self.assertTrue(solver.stats()["t_wall_total"]>=0)

    integr_options = {}
    integr_options["simplify"] = True
    integr_options["simplify_options"] = opts
    solver = ca.integrator("integrator","rk",{"x":x,"ode":f},integr_options)
    solver()
    self.assertTrue("t_proc_total" in solver.stats())
    self.assertTrue(solver.stats()["t_wall_total"]>=0)

    A = ca.DM.rand(3,3)
    if ca.has_linsol("lapacklu"):
      solver = ca.Linsol("linsol","lapacklu",A.sparsity(),opts)
      solver.solve(A,ca.vertcat(1,2,3))
      self.assertTrue("t_proc_total" in solver.stats())
      self.assertTrue(solver.stats()["t_wall_total"]>=0)

  def test_unicode(self):
    import sys
    import shutil
    if sys.version_info[0] < 3: return
    path = "üni/tüst.casadi"
    path2 = "üni/stüst.casadi"
    if "ghc-filesystem" not in ca.CasadiMeta.feature_list():
        path = "tüst.casadi"
        path2 = "stüst.casadi"
    x = ca.MX.sym("x")
    f = ca.Function("f",[x],[x**2])
    f.save(path)
    shutil.copy(path,path2)
    
    f2 = ca.Function.load(path2)
    self.checkfunction_light(f,f2,inputs=[3])
    
    
if __name__ == '__main__':
    unittest.main()
