from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *

scipy_available = True
try:
    import scipy.special
    from scipy.linalg import expm
except:
    scipy_available = False

class Misctests(casadiTestCase):
    
  def test_issue179A(self):
    self.message('Regression test #179 (A)')
    x = SX("x")
    f = SXFunction([x], [2 * x])
    f.init()
    y = f.eval([x])[0].data()
    print y
    
  def test_issue179B(self):
    self.message('Regression test #179 (B)')
    def calc_sparsity():
      x = casadi.SX("x")
      f = casadi.SXFunction([[x]], [[x ** 2]])
      f.init()
      return f.jacSparsity()
    
    def print_sparsity():
        sparsity = calc_sparsity()
        str(sparsity) # Segfault
        
    print_sparsity()
    
  def test_sanity(self):
    DMatrix(3,4,[1,2,1],[0,2,2,3],[0.738,0.39,0.99])
    self.assertRaises(RuntimeError,lambda : DMatrix(4,4,[1,2,1],[0,2,2,3],[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DMatrix(3,4,[1,2,1],[0,2,2,12],[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DMatrix(3,4,[1,2,1],[-10,2,2,3],[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DMatrix(3,4,[8,2,1],[0,2,2,3],[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DMatrix(3,4,[-3,2,1],[0,2,2,3],[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DMatrix(3,4,[1,2,1,2],[0,2,2,3],[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DMatrix(3,4,[1,2,1],[0,2,0,3],[0.738,0.39,0.99]))
  
  def test_setoptionerrors(self):
    self.message("option errors")
    x = SX("x")
    f = SXFunction([x],[x])
    
    f.setOption("name","foobar")
    self.assertRaises(RuntimeError,lambda : f.getOption("foobar"))
    self.assertRaises(RuntimeError,lambda : f.setOption("foobar",123))
    self.assertRaises(RuntimeError,lambda : f.setOption("name",123))
    
    self.assertRaises(RuntimeError,lambda : f.setOption("ad_mode","foo"))
    
    x = SX("x")
    f = SXFunction([x],[x])
    
    print "IpoptSolver"
    g = IpoptSolver(f)
    
    self.assertRaises(RuntimeError,lambda : g.setOption("monitor",["abc"]))
    g.setOption("monitor",["eval_f"])
    
    
  def test_copyconstr_norefcount(self):
    self.message("Copy constructor for non-refcounted classes")
    x = DMatrix(2,3,1)

    y = DMatrix(x)
    x[0,0] = 5
    
    self.assertFalse(id(x)==id(y))
    self.assertEqual(x[0,0],5)
    self.assertEqual(y[0,0],1)
    
  def test_copyconstr_refcount(self):
    self.message("Copy constructor for refcounted classes")
    x = sp_diag(4)

    y = CRSSparsity(x)
        
    x.resize(2,8)
    
    self.assertFalse(id(x)==id(y))
    self.assertTrue(x.numel(),y.numel())
    self.checkarray(x.shape,(2,8),"shape")
    self.checkarray(y.shape,(4,4),"shape")
    

  def test_copyconstr_refcount_lazy(self):
    self.message("Copy constructor for refcounted classes - lazy")
    x = SX("x")

    f = SXFunction([x],[2*x])
    f.init()
    f.input(0).setAll(2)
    g = SXFunction(f)

    f.input(0).set(5)
    f.evaluate()

    self.assertEqual(g.input(0),5)
    self.assertEqual(g.output(),10)

    
  def test_copy_norefcount(self):
    self.message("Shallow copy for non-refcounted classes")
    import copy
    
    x = DMatrix(2,3,1)

    y = copy.copy(x)
    x[0,0] = 5
    
    self.assertFalse(id(x)==id(y))
    self.assertEqual(x[0,0],5)
    self.assertEqual(y[0,0],1)
    
  def test_copy_refcount(self):
    self.message("Shallow copy for refcounted classes")
    import copy
    x = sp_diag(4)

    y = copy.copy(x)
        
    x.resize(2,8)
    
    self.assertFalse(id(x)==id(y))
    self.assertTrue(x.numel(),y.numel())
    self.checkarray(x.shape,(2,8),"shape")
    self.checkarray(y.shape,(4,4),"shape")
    
  def test_copy_refcount_lazy(self):
    self.message("Shallow copy for refcounted classes - lazy")
    import copy
    x = SX("x")

    f = SXFunction([x],[2*x])
    f.init()
    f.input(0).setAll(2)
    g = copy.copy(f)

    f.input(0).set(5)
    f.evaluate()

    self.assertEqual(g.input(0),5)
    self.assertEqual(g.output(),10)
    
  def test_deepcopy_norefcount(self):
    self.message("Deep copy for non-refcounted classes")
    import copy
    
    x = DMatrix(2,3,1)

    y = copy.deepcopy(x)
    x[0,0] = 5
    
    self.assertFalse(id(x)==id(y))
    self.assertEqual(x[0,0],5)
    self.assertEqual(y[0,0],1)
    
  def test_deepcopy_refcount(self):
    self.message("Deep copy for refcounted classes")
    import copy
    x = sp_diag(4)

    y = copy.deepcopy(x)
        
    x.resize(2,8)
    
    self.assertFalse(id(x)==id(y))
    self.assertTrue(x.numel(),y.numel())
    self.checkarray(x.shape,(2,8),"shape")
    self.checkarray(y.shape,(4,4),"shape")
    
  def test_deepcopy_refcount_lazy(self):
    self.message("Deep copy for refcounted classes - lazy")
    import copy
    x = SX("x")

    f = SXFunction([x],[2*x])
    f.init()
    f.input(0).setAll(2)
    g = copy.deepcopy(f)

    f.input(0).set(5)
    f.evaluate()

    self.assertEqual(g.input(0),2)
    self.assertEqual(g.output(),0)
    
  def test_options_introspection(self):
    self.message("options introspection")
    x=SX("x")
    f = SXFunction([x],[x**2])
    f.init()
    i = IpoptSolver(f)
    
    opts = i.getOptionNames()
    self.assertTrue(isinstance(opts,tuple))
    
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
    
if __name__ == '__main__':
    unittest.main()
    
