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
import pickle

scipy_available = True
try:
    import scipy.special
    from scipy.linalg import expm
except:
    scipy_available = False

class Misctests(casadiTestCase):
    
  def test_issue179A(self):
    self.message('Regression test #179 (A)')
    x = SXElement.sym("x")
    f = SXFunction([x], [2 * x])
    f.init()
    y = f.call([x])[0].data()
    print y
    
  def test_issue179B(self):
    self.message('Regression test #179 (B)')
    def calc_sparsity():
      x = casadi.SXElement.sym("x")
      f = casadi.SXFunction([x], [x ** 2])
      f.init()
      return f.jacSparsity()
    
    def print_sparsity():
        sparsity = calc_sparsity()
        str(sparsity) # Segfault
        
    print_sparsity()
    
  def test_sanity(self):
    DMatrix(Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.39,0.99])
    self.assertRaises(RuntimeError,lambda : DMatrix(Sparsity(4,4,[0,2,2,3],[1,2,1]),[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DMatrix(Sparsity(4,3,[0,2,2,12],[1,2,1]),[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DMatrix(Sparsity(4,3,[-10,2,2,3],[1,2,1]),[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DMatrix(Sparsity(4,3,[0,2,2,3],[8,2,1]),[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DMatrix(Sparsity(4,3,[0,2,2,3],[-3,2,1]),[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DMatrix(Sparsity(4,3,[0,2,2,3],[1,2,1,2]),[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DMatrix(Sparsity(4,3,[0,2,0,3],[1,2,1]),[0.738,0.39,0.99]))
  
  def test_setoptionerrors(self):
    self.message("option errors")
    x = SXElement.sym("x")
    f = SXFunction([x],[x])
    
    f.setOption("name","foobar")
    self.assertRaises(RuntimeError,lambda : f.getOption("foobar"))
    self.assertRaises(RuntimeError,lambda : f.setOption("foobar",123))
    self.assertRaises(RuntimeError,lambda : f.setOption("name",123))
    
    self.assertRaises(RuntimeError,lambda : f.setOption("ad_mode","foo"))
    
    x = SX.sym("x")
    nlp = SXFunction(nlpIn(x=x),nlpOut(f=x))

    try:
        print "ipopt"
        g = NlpSolver("ipopt", nlp)
    except:
        return
    
    self.assertRaises(RuntimeError,lambda : g.setOption("monitor",["abc"]))
    self.assertRaises(RuntimeError,lambda : g.setOption("monitor",["eval_f","abc"]))
    g.setOption("monitor",["eval_f"])
    
    
  def test_copyconstr_norefcount(self):
    self.message("Copy constructor for non-refcounted classes")
    x = DMatrix.ones(2,3)

    y = DMatrix(x)
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
    

  def test_copyconstr_refcount_lazy(self):
    self.message("Copy constructor for refcounted classes - lazy")
    x = SXElement.sym("x")

    f = SXFunction([x],[2*x])
    f.init()
    f.setInput(2,0)
    g = SXFunction(f)

    f.setInput(5,0)
    f.evaluate()

    self.assertEqual(g.getInput(0),5)
    self.assertEqual(g.getOutput(),10)

    
  def test_copy_norefcount(self):
    self.message("Shallow copy for non-refcounted classes")
    import copy
    
    x = DMatrix.ones(2,3)

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
    
  def test_copy_refcount_lazy(self):
    self.message("Shallow copy for refcounted classes - lazy")
    import copy
    x = SXElement.sym("x")

    f = SXFunction([x],[2*x])
    f.init()
    f.setInput(2,0)
    g = copy.copy(f)

    f.setInput(5,0)
    f.evaluate()

    self.assertEqual(g.getInput(0),5)
    self.assertEqual(g.getOutput(),10)
    
  def test_deepcopy_norefcount(self):
    self.message("Deep copy for non-refcounted classes")
    import copy
    
    x = DMatrix.ones(2,3)

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
    
  def test_deepcopy_refcount_lazy(self):
    self.message("Deep copy for refcounted classes - lazy")
    import copy
    x = SXElement.sym("x")

    f = SXFunction([x],[2*x])
    f.init()
    f.setInput(2,0)
    g = copy.deepcopy(f)

    f.setInput(5,0)
    f.evaluate()

    self.assertEqual(g.getInput(0),2)
    self.assertEqual(g.getOutput(),0)

  @requiresPlugin(NlpSolver,"ipopt")
  def test_options_introspection(self):
    self.message("options introspection")
    x=SX.sym("x")
    nlp = SXFunction(nlpIn(x=x),nlpOut(f=x**2))
    i = NlpSolver("ipopt", nlp)
    
    opts = i.getOptionNames()
    self.assertTrue(isinstance(opts,list))
    
    n = opts[0]
    self.assertTrue(type(n)==type(""))
    
    n = "monitor"
    
    d = i.getOptionDescription(n)
    self.assertTrue(type(d)==type(""))
    self.assertTrue(not("d"=="N/A"))
    
    d = i.getOptionTypeName(n)
    self.assertEqual(d,"OT_STRINGVECTOR")

    #d = i.getOptionAllowed(n)
    
  def test_monotonicity(self):
    self.message("monotonicity tests")
    l = []
    self.assertTrue(isIncreasing(l))
    self.assertTrue(isDecreasing(l))
    self.assertTrue(isNonIncreasing(l))
    self.assertTrue(isNonDecreasing(l))
    self.assertTrue(isMonotone(l))
    self.assertTrue(isStrictlyMonotone(l))
    l = [-3]
    self.assertTrue(isIncreasing(l))
    self.assertTrue(isDecreasing(l))
    self.assertTrue(isNonIncreasing(l))
    self.assertTrue(isNonDecreasing(l))
    self.assertTrue(isMonotone(l))
    self.assertTrue(isStrictlyMonotone(l))
    for l in [ [3,5], [-Inf,5], [3,Inf] , [-Inf,Inf] ]:
      self.assertTrue(isIncreasing(l))
      self.assertFalse(isDecreasing(l))
      self.assertFalse(isNonIncreasing(l))
      self.assertTrue(isNonDecreasing(l))
      self.assertTrue(isMonotone(l))
      self.assertTrue(isStrictlyMonotone(l))
    l = [5,3]
    self.assertFalse(isIncreasing(l))
    self.assertTrue(isDecreasing(l))
    self.assertTrue(isNonIncreasing(l))
    self.assertFalse(isNonDecreasing(l))
    self.assertTrue(isMonotone(l))
    self.assertTrue(isStrictlyMonotone(l))
    for l in [ [5,5], [-Inf,-Inf], [Inf,Inf] ]:
      self.assertFalse(isIncreasing(l))
      self.assertFalse(isDecreasing(l))
      self.assertTrue(isNonIncreasing(l))
      self.assertTrue(isNonDecreasing(l))
      self.assertTrue(isMonotone(l))
      self.assertFalse(isStrictlyMonotone(l))
    for l in [ [3,2,5], [1,NaN,2], [NaN] ] :
      self.assertFalse(isIncreasing(l))
      self.assertFalse(isDecreasing(l))
      self.assertFalse(isNonIncreasing(l))
      self.assertFalse(isNonDecreasing(l))
      self.assertFalse(isMonotone(l))
      self.assertFalse(isStrictlyMonotone(l))
    
  def test_regression418(self):
    self.message("Segfault regression check")
    f = ControlSimulator()
    try:
      f.setOption("integrator_options",None) # This should not give a segfault
    except:
      pass
  
  def test_regression448(self):
    self.message("regression test for segfaukt when printing")
    x = SX.sym("x")

    f = SXFunction(controldaeIn(x=x),daeOut(ode=x))
    f.init()

    sim = ControlSimulator(f,[0,1])
    sim.setOption("integrator_options",{"abstol": 1e-4})
    sim.printOptions()
    
  def test_IOscheme_indexing(self):
    self.message("IOscheme indexing")
    s = daeIn(x=2)
    
    self.assertEqual(s[0],2)
    self.assertEqual(s["x"],2)
    with self.assertRaises(Exception):
      s["xfgfd"]
    with self.assertRaises(Exception):
      s[100]
    s[-1]
    
  def test_IOscheme_output(self):
    x = 2
    p = 3
    s = daeIn(x=x,p=p)
    
    pp, = daeIn(s,"p")
    self.assertEqual(pp,p)
    
    xx,pp = daeIn(s,"x","p")
    self.assertEqual(pp,p)
    self.assertEqual(xx,x)
    
  def test_pickling(self):

    a = Sparsity.tril(4)
    s = pickle.dumps(a)
    b = pickle.loads(s)
    self.assertTrue(a==b)

    a = Sparsity()
    s = pickle.dumps(a)
    b = pickle.loads(s)
    self.assertTrue(a.isNull())
    
    a = IMatrix(Sparsity.tril(4),range(10))
    s = pickle.dumps(a)
    b = pickle.loads(s)
    self.checkarray(a,b)


    a = DMatrix(Sparsity.tril(4),range(10))
    s = pickle.dumps(a)
    b = pickle.loads(s)
    self.checkarray(a,b)
    
  def test_exceptions(self):
    try:
      NlpSolver(123)
      self.assertTrue(False)
    except NotImplementedError as e:
      print e.message
      assert "NlpSolver(str,Function)" in e.message
      assert "You have: NlpSolver(int)" in e.message
      assert "::" not in e.message
      assert "std" not in e.message

    try:
      vertcat(123)
      self.assertTrue(False)
    except NotImplementedError as e:
      print e.message
      assert "vertcat([SX]" in e.message
      assert "vertcat([array(double)" in e.message
      assert "You have: vertcat(int)" in e.message
      assert "::" not in e.message
      assert "std" not in e.message
      
    try:
      substitute(123)
      self.assertTrue(False)
    except NotImplementedError as e:
      print e.message
      assert "substitute(SX,SX,SX)" in e.message
      assert "substitute([SX] ,[SX] ,[SX] )" in e.message
      assert "You have: substitute(int)" in e.message
      assert "::" not in e.message
      assert "std" not in e.message
      
    try:
      SXFunction(daeIn(x=SX.sym("x")))
      self.assertTrue(False)
    except NotImplementedError as e:
      print e.message
      assert "SXFunction(scheme(SX),[SX] )" in e.message
      assert "You have: SXFunction(scheme(SX))" in e.message
      assert "::" not in e.message
      assert "std" not in e.message

    try:
      NlpSolver.loadPlugin(132)
      self.assertTrue(False)
    except TypeError as e:
      print e.message
      assert "type 'str' expected" in e.message
      assert "NlpSolver.loadPlugin" in e.message
      assert "You have: NlpSolver.loadPlugin(int)" in e.message
      assert "::" not in e.message
      assert "std" not in e.message

    x=SX.sym("x")
      
    try:
      [x]+ x
      self.assertTrue(False)
    except TypeError as e:
      print e.message
      assert "You try to do: [SX] + SX" in e.message

    try:
      x + [x]
      self.assertTrue(False)
    except TypeError as e:
      print e.message
      assert "You try to do: SX + [SX]" in e.message


    try:
      daeIn(x=x,p=[x])
      self.assertTrue(False)
    except TypeError as e:
      print e.message
      assert "You have: (x=SX, p=[SX])" in e.message

    try:
      QpSolver("qp",123)
      self.assertTrue(False)
    except NotImplementedError as e:
      print e.message
      assert "QpSolver(str,QPStructure)" in e.message
      
    try:
      SXFunction(qpStruct(a=12),[x])
      self.assertTrue(False)
    except NotImplementedError as e:
      print e.message
      assert "QPStructure([Sparsity] )" in e.message
      assert "You have: QPStructure([int,Sparsity])" in e.message
      assert "QPStructure(a=int)" in e.message
  
    try:
      x.reshape(2)
      self.assertTrue(False)
    except NotImplementedError as e:
      print e.message
      assert "reshape(SX,(int,int) )" in e.message

    try:
      x.reshape(("a",2))
      self.assertTrue(False)
    except NotImplementedError as e:
      print e.message
      assert "You have: reshape((str,int))" in e.message
      
    try:
      diagsplit("s")
      self.assertTrue(False)
    except NotImplementedError as e:
      print e.message
      assert "diagsplit(SX ,int)" in e.message
      assert "diagsplit(array(double) ,int)" in e.message
      
  def test_callkw(self):
      x = SX.sym("x")

      f = SXFunction(nlpIn(x=x),nlpOut(g=x**2))
      f.init()

      [f_,g_] = f(x=4)
      self.checkarray(g_,DMatrix(16))

      with self.assertRaises(RuntimeError):
        [f_,g_] = f(m=4)
      
      try:
        [f_,g_] = f(x=Sparsity.dense(2))
        self.assertTrue(False)
      except RuntimeError as e:
        self.assertTrue("Function(scheme(SX))" in e.message)
        self.assertTrue("Function([SX] )" in e.message)
        self.assertTrue("You have: Function(scheme(Sparsity))" in e.message)

      with self.assertRaises(RuntimeError):
        [f_,g_] = f(x=[x])


      f = SXFunction([x],nlpOut(g=x**2))
      f.init()

      with self.assertRaises(Exception):
        [f_,g_] = f(x=4)
        
  def test_getscheme(self):
    x = SX.sym("x")
    p = SX.sym("p")

    F = SXFunction(nlpIn(x=x,p=p),nlpOut(g=x**2,f=x+p))
    F.init()
    
    fc = F(x=3,p=4)
    [f] = fc.get.f
    self.checkarray(f,DMatrix([7]))
    [g] = fc.get.g
    self.checkarray(g,DMatrix([9]))
    [f,g] = fc.get.f.g
    self.checkarray(f,DMatrix([7]))
    self.checkarray(g,DMatrix([9]))
    [g,f] = fc.get.g.f
    self.checkarray(f,DMatrix([7]))
    self.checkarray(g,DMatrix([9]))
    
  def test_assertions(self):
    
    x = MX.sym("x") 
    
    z = x**2
    
    z = z.attachAssert(z>3,"x must be larger than 3")
    
    v = sin(z)
    
    f = MXFunction([x],[v])
    f.init()
    
    print f
    
    f.setInput(-6)
    f.evaluate()
    
    f.setInput(1)
    
    try :
      f.evaluate()
    except Exception as e:
      print str(e)
      self.assertTrue("x must be larger than 3" in str(e))
    

    
pickle.dump(Sparsity(),file("temp.txt","w"))
    
if __name__ == '__main__':
    unittest.main()
    
