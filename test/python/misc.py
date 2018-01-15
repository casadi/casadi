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
import pickle
from operator import itemgetter
import sys


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
      x = casadi.SX.sym("x")
      f = casadi.Function('f', [x], [x ** 2])
      return f.sparsity_jac(0, 0)

    def print_sparsity():
        sparsity = calc_sparsity()
        str(sparsity) # Segfault

    print_sparsity()

  def test_sanity(self):
    DM(Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.39,0.99])
    self.assertRaises(RuntimeError,lambda : DM(Sparsity(4,4,[0,2,2,3],[1,2,1]),[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DM(Sparsity(4,3,[0,2,2,12],[1,2,1]),[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DM(Sparsity(4,3,[-10,2,2,3],[1,2,1]),[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DM(Sparsity(4,3,[0,2,2,3],[8,2,1]),[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DM(Sparsity(4,3,[0,2,2,3],[-3,2,1]),[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DM(Sparsity(4,3,[0,2,2,3],[1,2,1,2]),[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DM(Sparsity(4,3,[0,2,0,3],[1,2,1]),[0.738,0.39,0.99]))

  def test_copyconstr_norefcount(self):
    self.message("Copy constructor for non-refcounted classes")
    x = DM.ones(2,3)

    y = DM(x)
    x[0,0] = 5

    self.assertFalse(id(x)==id(y))
    self.assertEqual(x[0,0],5)
    self.assertEqual(y[0,0],1)

  def test_copyconstr_refcount(self):
    self.message("Copy constructor for refcounted classes")
    x = Sparsity.diag(4)

    y = Sparsity(x)

    x.resize(2,8)

    self.assertFalse(id(x)==id(y))
    self.assertTrue(x.numel(),y.numel())
    self.checkarray(x.shape,(2,8),"shape")
    self.checkarray(y.shape,(4,4),"shape")

  def test_copy_norefcount(self):
    self.message("Shallow copy for non-refcounted classes")
    import copy

    x = DM.ones(2,3)

    y = copy.copy(x)
    x[0,0] = 5

    self.assertFalse(id(x)==id(y))
    self.assertEqual(x[0,0],5)
    self.assertEqual(y[0,0],1)

  def test_copy_refcount(self):
    self.message("Shallow copy for refcounted classes")
    import copy
    x = Sparsity.diag(4)

    y = copy.copy(x)

    x.resize(2,8)

    self.assertFalse(id(x)==id(y))
    self.assertTrue(x.numel(),y.numel())
    self.checkarray(x.shape,(2,8),"shape")
    self.checkarray(y.shape,(4,4),"shape")

  def test_deepcopy_norefcount(self):
    self.message("Deep copy for non-refcounted classes")
    import copy

    x = DM.ones(2,3)

    y = copy.deepcopy(x)
    x[0,0] = 5

    self.assertFalse(id(x)==id(y))
    self.assertEqual(x[0,0],5)
    self.assertEqual(y[0,0],1)

  def test_deepcopy_refcount(self):
    self.message("Deep copy for refcounted classes")
    import copy
    x = Sparsity.diag(4)

    y = copy.deepcopy(x)

    x.resize(2,8)

    self.assertFalse(id(x)==id(y))
    self.assertTrue(x.numel(),y.numel())
    self.checkarray(x.shape,(2,8),"shape")
    self.checkarray(y.shape,(4,4),"shape")

  @requiresPlugin(nlpsol,"ipopt")
  def test_options_introspection(self):
    self.message("options introspection")
    x=SX.sym("x")
    nlp = {'x':x, 'f':x**2}
    i = nlpsol('i', "ipopt", nlp)

    opts = i.optionNames()
    self.assertTrue(isinstance(opts,list))

    n = opts[0]
    self.assertTrue(type(n)==type(""))

    n = "monitor"

    d = i.optionDescription(n)
    self.assertTrue(type(d)==type(""))
    self.assertTrue(not("d"=="N/A"))

    d = i.optionTypeName(n)
    self.assertEqual(d,"OT_STRINGVECTOR")

    #d = i.optionAllowed(n)

  @unittest.skipIf(sys.version_info>=(3,0),"too lazy to fix now")
  def test_pickling(self):

    a = Sparsity.lower(4)
    s = pickle.dumps(a)
    b = pickle.loads(s)
    self.assertTrue(a==b)

    a = Sparsity()
    s = pickle.dumps(a)
    b = pickle.loads(s)
    self.assertTrue(a.is_null())

    a = IM(Sparsity.lower(4),list(range(10)))
    s = pickle.dumps(a)
    b = pickle.loads(s)
    self.checkarray(a,b)


    a = DM(Sparsity.lower(4),list(range(10)))
    s = pickle.dumps(a)
    b = pickle.loads(s)
    self.checkarray(a,b)

  @known_bug()
  def test_exceptions(self):
    try:
      nlpsol(123)
      self.assertTrue(False)
    except NotImplementedError as e:
      print(e.message)
      assert "nlpsol(str,str,Function,Dict)" in e.message
      assert "You have: nlpsol(int)" in e.message
      assert "::" not in e.message
      assert "std" not in e.message

    try:
      vertcat(*123)
      self.assertTrue(False)
    except NotImplementedError as e:
      print(e.message)
      assert "vertcat(*[SX]" in e.message
      assert "vertcat(*[DM" in e.message
      assert "You have: vertcat(*int)" in e.message
      assert "::" not in e.message
      assert "std" not in e.message

    try:
      substitute(123)
      self.assertTrue(False)
    except NotImplementedError as e:
      print(e.message)
      assert "substitute(SX,SX,SX)" in e.message
      assert "substitute([SX] ,[SX] ,[SX] )" in e.message
      assert "You have: substitute(int)" in e.message
      assert "::" not in e.message
      assert "std" not in e.message

    try:
      load_nlpsol(132)
      self.assertTrue(False)
    except TypeError as e:
      print(e.message)
      assert "Failed to convert input to str" in e.message
      assert "::" not in e.message
      assert "std" not in e.message

    x=SX.sym("x")

    try:
      [x]+ x
      self.assertTrue(False)
    except TypeError as e:
      print(e.message)

    try:
      x + [x]
      self.assertTrue(False)
    except TypeError as e:
      print(e.message)

    try:
      x.reshape(2)
      self.assertTrue(False)
    except NotImplementedError as e:
      print(e.message)
      assert "reshape(SX,(int,int) )" in e.message

    try:
      x.reshape(("a",2))
      self.assertTrue(False)
    except NotImplementedError as e:
      print(e.message)
      assert "You have: reshape((str,int))" in e.message

    try:
      diagsplit("s")
      self.assertTrue(False)
    except NotImplementedError as e:
      print(e.message)
      assert "diagsplit(SX,int)" in e.message
      assert "diagsplit(DM ,int)" in e.message

    try:
      DM("df")
      self.assertTrue(False)
    except NotImplementedError as e:
      print(e.message)
      assert "  DM (" in e.message

    try:
      vertcat(*[1,SX.sym('x'),MX.sym('x')])
      self.assertTrue(False)
    except NotImplementedError as e:
      print(e.message)
      assert "  vertcat(*" in e.message

  def test_getscheme(self):
    x = SX.sym("x")
    p = SX.sym("p")

    F = Function('F', [x, p], [x+p, x**2], ['x', 'p'], ['f', 'g'])

    fc = F(x=3,p=4)
    f = fc['f']
    self.checkarray(f,DM([7]))
    g = fc['g']
    self.checkarray(g,DM([9]))
    [f,g] = itemgetter('f','g')(fc)
    self.checkarray(f,DM([7]))
    self.checkarray(g,DM([9]))
    [g,f] = itemgetter('g','f')(fc)
    self.checkarray(f,DM([7]))
    self.checkarray(g,DM([9]))

  def test_assertions(self):

    x = MX.sym("x")

    z = x**2

    z = z.attachAssert(z>3,"x must be larger than 3")

    v = sin(z)

    f = Function('f', [x],[v])

    print(f)

    f_out = f(-6)
    try :
      f_out = f(1)
    except Exception as e:
      print(str(e))
      self.assertTrue("x must be larger than 3" in str(e))

  @requires_nlpsol("ipopt")
  def test_output(self):
    with capture_stdout() as result:
      DM([1,2]).print_dense()


    assert "2" in result[0]

    x=SX.sym("x")
    f = {'x':x, 'f':x**2}
    solver = nlpsol("solver", "ipopt",f)
    with capture_stdout() as result:
      solver_out = solver(x0=0)

    assert "Number of nonzeros in equality constraint" in result[0]
    assert "iter    objective    inf_pr" in result[0]

    with capture_stdout() as result:
      try:
        solver = nlpsol("solver","foo",f)
      except:
        pass

    assert "casadi_nlpsol_foo" in result[1]


if __name__ == '__main__':
    unittest.main()
