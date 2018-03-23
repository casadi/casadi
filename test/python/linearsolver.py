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
import random

warnings.filterwarnings("ignore",category=DeprecationWarning)

lsolvers = []
try:
  load_linsol("csparse")
  lsolvers.append(("csparse",{},set()))
except:
  pass

"""
try:
  load_linsol("csparsecholesky")
  lsolvers.append(("csparsecholesky",{},{"posdef","symmetry"}))
except:
  pass
"""
try:
  load_linsol("lapacklu")
  lsolvers.append(("lapacklu",{},set()))
except:
  pass

try:
  load_linsol("lapackqr")
  lsolvers.append(("lapackqr",{},set()))
except:
  pass

try:
  load_linsol("symbolicqr")
  lsolvers.append(("symbolicqr",{},set()))
except:
  pass

try:
  load_linsol("lsqr")
  lsolvers.append(("lsqr",{},set()))
except:
  pass

try:
  load_linsol("ma27")
  lsolvers.append(("ma27",{},{"symmetry"}))
except:
  pass

try:
  load_linsol("qr")
  lsolvers.append(("qr",{},set()))
except:
  pass

try:
  load_linsol("ldl")
  lsolvers.append(("ldl",{},{"posdef","symmetry"}))
except:
  pass

nsolvers = []

def nullspacewrapper(name, sp, options):
  a = SX.sym("a",sp)
  f = Function(name, [a],[nullspace(a)],options)
  return f

nsolvers.append((nullspacewrapper,{},set()))

print("linear solvers", lsolvers)

class LinearSolverTests(casadiTestCase):

  def test_nullspace(self):

    for A in  [
                  DM([[1,1.3],[2,5],[1,0.5],[1.8,1.7]]),
                  DM([[1,1.3],[2,5],[1,0.5]]),
                  DM([[1,1.3],[2,5],[1,0.5],[0.2,0.3],[-0.3,0.7]]),
                  DM([[1,0],[0,0],[0,1],[0,0]]),
                  DM([[1.3,0,0.4,1],[0.2,0.1,11,0],[0,1,0,0],[0.7,0.9,0,0],[1.1,0.99,0,0]])
              ]:
      n ,m = A.shape
      for Solver, options,req in nsolvers:
        if "symmetry" in req: continue
        solver = Solver("solver", A.T.sparsity(), options)
        solver_in = [0]*solver.n_in();solver_in[0]=A.T

        solver_out = solver.call(solver_in)

        self.checkarray(mtimes(A.T,solver_out[0]),DM.zeros(m,n-m))
        self.checkarray(mtimes(solver_out[0].T,solver_out[0]),DM.eye(n-m))

        options["ad_weight"] = 0
        options["ad_weight_sp"] = 0
        solver = Solver("solver", A.T.sparsity(), options)

        Jf = solver.jacobian_old(0, 0)

        options["ad_weight"] = 1
        options["ad_weight_sp"] = 1
        solver = Solver("solver", A.T.sparsity(), options)

        Jb = solver.jacobian_old(0, 0)

        Jf_in = [0]*Jf.n_in();Jf_in[0]=A.T
        Jb_in = [0]*Jb.n_in();Jb_in[0]=A.T

        Jf_out = Jf.call(Jf_in)
        Jb_out = Jb.call(Jb_in)

        self.checkarray(Jf_out[0],Jb_out[0])
        self.checkarray(Jf_out[1],Jb_out[1])

        d = solver.forward(1)

        r = numpy.random.rand(*A.shape)

        d_in = [A.T, 5, r.T]

        d_out = d.call(d_in)

        exact = d_out[0]

        solver_in = [0]*solver.n_in();solver_in[0]=A.T
        solver_out = solver(*solver_in)
        nom = solver_out

        eps = 1e-6
        solver_in = [0]*solver.n_in();solver_in[0]=(A+eps*r).T
        solver_out = solver(*solver_in)
        pert = solver_out

        fd = (pert-nom)/eps

        #print exact, fd

        #print numpy.linalg.svd(horzcat(*[exact, fd]).T)[1]

        #print "fd:", mtimes(fd.T,fd), numpy.linalg.eig(mtimes(fd.T,fd))[0]
        #print "exact:", mtimes(exact.T,exact), numpy.linalg.eig(mtimes(exact.T,exact))[0]
        #print "fd:", mtimes(fd,fd.T), numpy.linalg.eig(mtimes(fd,fd.T))[0]
        #print "exact:", mtimes(exact,exact.T), numpy.linalg.eig(mtimes(exact,exact.T))[0]

        V = numpy.random.rand(A.shape[0]-A.shape[1],A.shape[0]-A.shape[1])
        V = V+V.T
        print(V)
        #V = DM.eye(A.shape[0]-A.shape[1])
        a = mtimes([nom,V,fd.T])+mtimes([fd,V,nom.T])
        b = mtimes([nom,V,exact.T])+mtimes([exact,V,nom.T])

        print("here:", a-b)

        #self.checkarray(a,b,digits=5)

        V = numpy.random.rand(A.shape[0],A.shape[0])
        V = V+V.T
        V = DM.eye(A.shape[0])
        a = mtimes([nom.T,V,fd])+mtimes([fd.T,V,nom])
        b = mtimes([nom.T,V,exact])+mtimes([exact.T,V,nom])

        self.checkarray(a,b,digits=5)

  def test_simple_solve(self):
    A_ = DM([[3,7],[1,2]])
    b_ = DM([1,0.5])

    A = MX.sym("A",A_.sparsity())
    b = MX.sym("b",b_.sparsity())

    for Solver, options,req in lsolvers:
      C = solve(A,b,Solver,options)
      if "symmetry" in req:
        A0 = A_.T+A_
      else:
        A0 = A_
      f = Function("f", [A,b],[C])
      f_out = f(A0, b_)

      self.checkarray(f_out,np.linalg.solve(A0,b_))
      self.checkarray(mtimes(A0,f_out),b_)

  def test_pseudo_inverse(self):
    numpy.random.seed(0)
    A_ = DM(numpy.random.rand(4,6))

    A = MX.sym("A",A_.sparsity())
    As = SX.sym("A",A_.sparsity())

    for Solver, options,req in lsolvers:
      if "symmetry" in req: continue
      print(Solver)
      B = pinv(A,Solver,options)

      f = Function("f", [A],[B])
      f_out = f(A_)

      self.checkarray(mtimes(A_,f_out),DM.eye(4))

      f = Function("f", [As],[pinv(As)])
      f_out = f(A_)

      self.checkarray(mtimes(A_,f_out),DM.eye(4))

      solve(mtimes(A,A.T),A,Solver,options)
      pinv(A_,Solver,options)

      #self.checkarray(mtimes(A_,pinv(A_,Solver,options)),DM.eye(4))

    A_ = DM(numpy.random.rand(3,5))

    A = MX.sym("A",A_.sparsity())
    As = SX.sym("A",A_.sparsity())

    for Solver, options,req in lsolvers:
      print(Solver)
      if "symmetry" in req: continue
      B = pinv(A,Solver,options)

      f = Function("f", [A],[B])
      f_out = f(A_)
      self.checkarray(mtimes(A_,f_out),DM.eye(3))

      f = Function("f", [As],[pinv(As)])
      f_out = f(A_)

      self.checkarray(mtimes(A_,f_out),DM.eye(3))

      #self.checkarray(mtimes(pinv(A_,Solver,options),A_),DM.eye(3))

  def test_simple_solve_dmatrix(self):
    A = DM([[3,7],[1,2]])
    b = DM([1,0.5])
    for Solver, options, req in lsolvers:
      print(Solver)
      if "symmetry" in req:
        A0 = A.T+A
      else:
        A0 = A
      C = solve(A0,b,Solver,options)

      sol = np.linalg.solve(A0,b)
      self.checkarray(C,sol)
      self.checkarray(mtimes(A0,sol),b)

  def test_simple_trans(self):
    A = DM([[3,1],[7,2]])
    for Solver, options,req in lsolvers:
      if "symmetry" in req:
        A0 = A.T+A
      else:
        A0 = A
      solver = casadi.Linsol("solver", Solver, A0.sparsity(), options)
      b = DM([1,0.5])
      x = solver.solve(A0.T, b)
      res = np.linalg.solve(A0.T,b)
      self.checkarray(x, res)

  def test_simple(self):
    A = DM([[3,1],[7,2]])
    for Solver, options, req in lsolvers:
      print(Solver)
      if "symmetry" in req:
        A0 = A.T+A
      else:
        A0 = A
      solver = casadi.Linsol("solver", Solver, A0.sparsity(), options)
      b = DM([1,0.5])
      x = solver.solve(A0, b)
      res = np.linalg.solve(A0,b)
      self.checkarray(x, res)

  def test_simple_function_indirect(self):

    for Solver, options,req in lsolvers:
      A_ = DM([[3,1],[7,2]])
      A = MX.sym("A",A_.sparsity())
      b_ = DM([1,0.5])
      b = MX.sym("b",b_.sparsity())
      if "symmetry" in req: A_ = A_.T + A_
      print(Solver)
      solver = casadi.Linsol("solver", Solver, A.sparsity(), options)
      solver_in = {}
      solver_in["A"]=A_
      solver_in["B"]=b_

      sol = {"X": solver.solve(A,b)}
      sol["A"] = A
      sol["B"] = b
      relay = Function("relay", sol, ["A","B"], ["X"])

      A_0 = A[0,0]
      A_1 = A[0,1]
      A_2 = A[1,0]
      A_3 = A[1,1]

      b_0 = b[0]
      b_1 = b[1]

      solution = Function("solution", {"A":A, "B":b, "X":vertcat(*[(((A_3/((A_0*A_3)-(A_2*A_1)))*b_0)+(((-A_1)/((A_0*A_3)-(A_2*A_1)))*b_1)),((((-A_2)/((A_0*A_3)-(A_2*A_1)))*b_0)+((A_0/((A_0*A_3)-(A_2*A_1)))*b_1))])}, ["A", "B"], ["X"])


      self.checkfunction(relay,solution,inputs=solver_in)

      if Solver in ["qr","ldl"]:
        self.check_codegen(relay,inputs=solver_in)

  @memory_heavy()
  def test_simple_solve_node(self):

    for Solver, options, req in lsolvers:

      for A_,b_ in [
                       (DM([[3,1],[7,2]]),DM([[1,0.3],[0.5,0.7]])),
                       (sparsify(DM([[3,0],[7,2]])),DM([[1,0.3],[0.5,0.7]])),
                       (DM([[3,1],[7,2]]),sparsify(DM([[1,0],[0,0.7]])))
                   ]:
        if "symmetry" in req: A_ = A_.T+A_

        A = MX.sym("A",A_.sparsity())
        b = MX.sym("b",b_.sparsity())
        print(Solver)
        solver = casadi.Linsol("solver", Solver, A.sparsity(), options)
        for tr in [True, False]:
          x = solver.solve(A,b,tr)
          f = Function("f", [A,b],[x])
          f_out = f(A_, b_)

          if tr:
            A_0 = A[0,0]
            A_1 = A[1,0]
            A_2 = A[0,1]
            A_3 = A[1,1]
          else:
            A_0 = A[0,0]
            A_1 = A[0,1]
            A_2 = A[1,0]
            A_3 = A[1,1]

          b_0 = b[0,0]
          b_1 = b[1,0]

          c_0 = b[0,1]
          c_1 = b[1,1]

          solution = Function("solution", [A,b],[blockcat([[(((A_3/((A_0*A_3)-(A_2*A_1)))*b_0)+(((-A_1)/((A_0*A_3)-(A_2*A_1)))*b_1)),(((A_3/((A_0*A_3)-(A_2*A_1)))*c_0)+(((-A_1)/((A_0*A_3)-(A_2*A_1)))*c_1))],[((((-A_2)/((A_0*A_3)-(A_2*A_1)))*b_0)+((A_0/((A_0*A_3)-(A_2*A_1)))*b_1)),((((-A_2)/((A_0*A_3)-(A_2*A_1)))*c_0)+((A_0/((A_0*A_3)-(A_2*A_1)))*c_1))]])])

          solution_in = [0]*solution.n_in();solution_in[0]=A_
          solution_in[1]=b_

          self.checkfunction(f,solution,inputs=solution_in)

          if "SymbolicQR" not in str(Solver) : continue
          solversx = f.expand('expand_'+f.name())
          solversx_in = [0]*solversx.n_in();solversx_in[0]=A_
          solversx_in[1]=b_

          self.checkfunction(solversx,solution,digits_sens = 7)

  def test_large_sparse(self):

    n = 10

    for Solver, options,req in lsolvers:
      numpy.random.seed(1)
      A = self.randDM(n,n,sparsity=0.5)
      if "symmetry" in req: A = A.T+A
      b = self.randDM(n,3,sparsity=0.5)

      As = MX.sym("A",A.sparsity())
      bs = MX.sym("B",b.sparsity())

      C = solve(A,b,Solver,options)

      self.checkarray(mtimes(A,C),b)

      f = Function("f", [As,bs],[solve(As,bs,Solver,options)])
      f_out = f(A, b)

      self.checkarray(mtimes(A,f_out),b)

  def test_ma27(self):
      n = np.nan

      A = DM([[2, 3, n, n, n],
             [3, 0, 4, n, 6],
             [n, 4, 1, 5, n],
             [n, n, 5, 0, n],
             [n, 6, n, n, 1]])
      A = A[sparsify(A==A).sparsity()]

      b = DM([8,45,31,15,17])

      ref = np.linalg.solve(A,b)
      for Solver, options, req in lsolvers:
        if "posdef" in req: continue
        As = MX.sym("A",A.sparsity())
        bs = MX.sym("B",b.sparsity())
        C = solve(A,b,Solver,options)
        self.checkarray(ref,C)


  def test_large_sparse2(self):
    numpy.random.seed(1)
    n = 10

    for Solver, options, req in lsolvers:

      print(Solver)
      A = self.randDM(n,n,sparsity=0.5)
      b = self.randDM(n,3)
      if "posdef" in req:
        A = mtimes(A.T, A)
        A = densify(A)
      elif "symmetry" in req:
        A = A.T+A
        A[Sparsity.diag(n)] =1e-8
        A = densify(A)

      As = MX.sym("A",A.sparsity())
      bs = MX.sym("B",b.sparsity())
      C = solve(A,b,Solver,options)
      digits = 7 if "ma" in str(Solver) else 10

      self.checkarray(mtimes(A,C),b,digits=digits)

      for As_,A_ in [(As,A),(densify(As),densify(A)),(densify(As).T,densify(A).T),(densify(As.T),densify(A.T)),(As.T,A.T)]:
        f = Function("f", [As,bs],[solve(As_,bs,Solver,options)])
        f_out = f(A, b)

        self.checkarray(mtimes(A_,f_out),b,digits=digits)

  def test_dimmismatch(self):
    A = DM.eye(5)
    b = DM.ones((4,1))
    As = MX.sym("A",A.sparsity())
    bs = MX.sym("b",b.sparsity())

    for Solver, options, req in lsolvers:
      with self.assertRaises(Exception):
        solve(A,b,Solver,options)

      with self.assertRaises(Exception):
        solve(As,bs,Solver,options)

  @memory_heavy()
  def test_thread_safety(self):
    x = MX.sym('x')
    y = MX.sym('y')


    for Solver, options, req in lsolvers:
      print(Solver)
      f = Function('f', [x, y], [x ** 2 - y])
      finv = rootfinder('finv', "newton", f, {"linear_solver": Solver, "linear_solver_options": options})

      finv_par = finv.map(200, 'thread',4)
      res = finv_par(numpy.ones(200), numpy.linspace(0, 10, 200))
      self.checkarray(norm_inf(res.T-sqrt(numpy.linspace(0, 10, 200))),0, digits=5)

      z = solve(blockcat([[1+x,0],[0,1+y]]),vertcat(x,y), Solver, options)
      f = Function('f',[x,y],[z])
      f_par = f.map(200, 'thread',4)
      res = f_par(numpy.linspace(10, 0, 200), numpy.linspace(0, 10, 200))




if __name__ == '__main__':
    unittest.main()
