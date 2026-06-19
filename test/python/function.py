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
import numpy as np
from numpy import inf, pi
import casadi as c
import numpy
import unittest
from types import *
from helpers import *
import pickle
import os
import re
import sys
from casadi.tools import capture_stdout
import glob
import gc
import tempfile

import random as random_base

scipy_interpolate = False
try:
  import scipy.interpolate
  scipy.interpolate.RectBivariateSpline
  scipy_interpolate = True
except:
  pass


class Functiontests(casadiTestCase):

  def test_call_empty(self):
    x = ca.SX.sym("x",2)
    fsx = ca.Function("fsx", [x,[]],[x])
    x = ca.MX.sym("x",2)
    fmx1 = ca.Function("fmx1", [x,ca.MX()],[x])
    fmx2 = ca.Function("fmx2", [x,[]],[x])

    for f in [fsx,fmx1,fmx2]:
      f(0,0)

      X = ca.MX.sym("X",2)
      F = f(X,ca.MX())
      g = ca.Function("g", [X],[F])

      g(0)

    x = ca.SX.sym("x",2)
    fsx = ca.Function("fsx", [x],[x,[]])
    x = ca.MX.sym("x",2)
    fmx1 = ca.Function("fmx1", [x],[x,ca.MX()])
    fmx2 = ca.Function("fmx2", [x],[x,[]])

    for f in [fsx,fmx1,]:
      f(0)

      X = ca.MX.sym("X",2)
      F = list(f(X))
      g = ca.Function("g", [X],F)

      g(0)

  def test_MX_funSeed(self):
    self.message("MX_funSeed")
    x1 = ca.MX.sym("x",2)
    y1 = ca.MX.sym("y")
    x2 = ca.MX.sym("x",2)
    y2 = ca.MX.sym("y")
    p= ca.Function("p", [x1,y1,x2,y2],[ca.sin(x1) + y1,ca.sin(x2) + y2])

    n1 = ca.DM([4,5])
    N1 = 3
    n2 = ca.DM([5,7])
    N2 = 8

    out = p(n1,N1,n2,N2)

    self.checkarray(ca.sin(n1)+N1,out[0],"output")
    self.checkarray(ca.sin(n2)+N2,out[1],"output")
  def test_segfault(self):
    f = ca.Function()
    with self.assertRaises(Exception):
      f.stats()
  def test_issue304(self):
    self.message("regression test for #304") # this code used to segfault
    x = ca.SX.sym("x")

    f = ca.Function("f", [x],[x**2,x**3])

    X = ca.MX.sym("X")

    z=list(f(X))

    g = ca.Function("g", [X], z).expand()

  def test_jacobian(self):
    x = ca.SX.sym("x",3,1)
    y = ca.SX.sym("y",2,1)

    f = ca.Function("f", [x,y],[x**2,y,x*y[0]])

    g = jacobian_old(f, 0, 0)

    self.assertEqual(g.n_in(),f.n_in())
    self.assertEqual(g.n_out(),f.n_out()+1)

  def test_xfunction(self):
    x = ca.SX.sym("x",3,1)
    y = ca.SX.sym("y",2,1)

    f = ca.Function("f", [x,y],[x**2,y,x*y[0]])

    X = ca.MX.sym("x",3,1)
    Y = ca.MX.sym("y",2,1)

    F = ca.Function("F", [X,Y],[X**2,Y,X*Y[0]])

    self.checkfunction(f,F,inputs=[[0.1,0.7,1.3],[7.1,2.9]],sens_der=False,evals=False)


  @memory_heavy()
  def test_jacobians(self):

    x = ca.SX.sym("x")

    self.assertEqual(ca.jacobian(5,x).nnz(),0)


    def test(sp):
      x = ca.SX.sym("x",sp.size2())
      sp2 = ca.jacobian(ca.DM.ones(sp) @ x,x).sparsity()
      self.checkarray(sp.row(),sp2.row());
      self.checkarray(sp.colind(),sp2.colind());

    for i in range(5):
      test(ca.Sparsity.lower(i))
      test(ca.Sparsity.lower(i).T)
      test(ca.Sparsity.dense(i,i))
      test(ca.Sparsity.diag(i))

    for i in [63,64,65,127,128,129]:
      d = ca.Sparsity.diag(i)
      test(d)

      test(d + ca.Sparsity.rowcol([0],[5],i,i))

      b = ca.Sparsity.band(i,-1) + ca.Sparsity.band(i,1)
      test(b + ca.Sparsity.rowcol([0],[5],i,i))

    m = ca.DM.ones(ca.Sparsity.diag(129))
    m[:50,0] = 1
    m[60:,0] = 1
    m[6:9,6] = 1
    m[9,9:12] = 1

    sp = m[:,:120].sparsity()

    test(sp)
    #test(sp.T)

    m = ca.DM.ones(ca.Sparsity.diag(64))
    m[:50,0] = 1
    m[60:,0] = 1

    sp = m.T[:40,:].sparsity()
    test(sp)
    test(sp.T)

    sp = m[:40,:].sparsity()
    test(sp)
    test(sp.T)

    sp = m.T[:20,:].sparsity()
    test(sp)
    test(sp.T)

    sp = m[:20,:].sparsity()
    test(sp)
    test(sp.T)

    for i in [63,64,65,127,128,129]:
      test(ca.Sparsity.lower(i))
      test(ca.Sparsity.lower(i).T)

    for n in ([63,64,65,127,128,129] if args.run_slow else [63,64,65]):
      for m in ([63,64,65,127,128,129] if args.run_slow else [63,64,65]):
        print((n,m))
        sp = ca.Sparsity.dense(n,m)

        test(sp)

        random.seed(0)

        I = ca.DM.ones(sp)
        for i in range(n):
          for j in range(m):
            if random.random()<0.5:
              I[i,j] = 0
        I = ca.sparsify(I)

        sp_holes = I.sparsity()

        test(sp_holes)

        z = ca.DM(sp_holes.size1(), sp_holes.size2())

        R = 5
        v = []
        for r in range(R):
          h = [z]*5
          h[r] = I
          v.append(ca.horzcat(*h))
        d = ca.vertcat(*v)

        test(d.sparsity())

  @memory_heavy()
  def test_hessians(self):
    def test(sp):
      x = ca.SX.sym("x",sp.size2())
      self.assertTrue(sp==sp.T)
      f = ca.Function("f", [x],[ca.mtimes([x.T,ca.DM.ones(sp),x])])
      J = hessian_old(f, 0, 0)
      sp2 = J.sparsity_out(0)
      self.checkarray(sp.row(),sp2.row())
      self.checkarray(sp.colind(),sp2.colind())

    A = ca.DM([[1,1,0,0,0,0],[1,1,1,0,1,1],[0,1,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,1,0,1,0,1]])
    A = ca.sparsify(A)
    C = A.sparsity()

    test(C)

    A = ca.DM([[1,0,0,0,0,0],[0,1,1,0,1,1],[0,1,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,1,0,1,0,1]])
    A = ca.sparsify(A)
    C = A.sparsity()

    test(C)

    A = ca.DM([[1,0,0,0,0,0],[0,1,0,0,1,1],[0,0,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,1,0,1,0,1]])
    A = ca.sparsify(A)
    C = A.sparsity()

    test(C)

    A = ca.DM([[0,0,0,0,0,0],[0,1,0,0,1,1],[0,0,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,1,0,1,0,1]])
    A = ca.sparsify(A)
    C = A.sparsity()

    test(C)

    A = ca.DM([[0,0,0,0,0,0],[0,1,0,0,1,0],[0,0,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,0,0,1,0,1]])
    A = ca.sparsify(A)
    C = A.sparsity()

    test(C)


    for i in [63,64,65,100,127,128,129]:
      d = ca.Sparsity.diag(i)
      test(d)

      test(d + ca.Sparsity.rowcol([0],[5],i,i) + ca.Sparsity.rowcol([5],[0],i,i))

      b = ca.Sparsity.band(i,-1) + ca.Sparsity.band(i,1)
      test(b)
      test(b + ca.Sparsity.rowcol([0],[5],i,i) + ca.Sparsity.rowcol([5],[0],i,i))

      d = ca.Sparsity.dense(i,i)
      test(d)

      d = ca.Sparsity.diag(i) + ca.Sparsity.triplet(i,i,[0]*i,list(range(i)))+ca.Sparsity.triplet(i,i,list(range(i)),[0]*i)
      test(d)


      sp = ca.Sparsity.dense(i,i)

      random.seed(0)

      I = ca.DM.ones(sp)
      for ii in range(i):
        for jj in range(i):
          if random.random()<0.5:
            I[ii,jj] = 0
            I[jj,ii] = 0
      I = ca.sparsify(I)

      sp_holes = I.sparsity()

      test(sp_holes)

      z = ca.DM(sp_holes.size1(), sp_holes.size2())

      R = 5
      v = []
      for r in range(R):
        h = [z]*5
        h[r] = I
        v.append(ca.horzcat(*h))
      d = ca.vertcat(*v)

      test(d.sparsity())

  def test_customIO(self):
    x = ca.SX.sym("x")
    f = ca.Function('f',[x],[x*x, x],["i0"], ["foo","bar"])

    ret = f(i0=12)

    self.checkarray(ca.DM([144]),ret["foo"])
    self.checkarray(ca.DM([12]),ret["bar"])


    with self.assertRaises(Exception):
      f_out["baz"]  # pyright: ignore[reportUndefinedVariable]

    ret = f(i0=ca.SX(12))
    self.checkarray(ret["foo"],ca.DM([144]))
    self.checkarray(ret["bar"],ca.DM([12]))
    with self.assertRaises(Exception):
      self.checkarray(ret["baz"],ca.DM([12]))

  def test_derivative_simplifications(self):

    n = 1
    x = ca.SX.sym("x",n)

    M = ca.Function("M", [x],[((x-ca.DM(list(range(n))))) @ x.T])

    P = ca.MX.sym("P",n,n)
    X = ca.MX.sym("X",n)

    M_X= M(X)

    Pf = ca.Function("P", [X, P], [M_X @ P])

    P_P = jacobian_old(Pf, 1, 0)

    self.assertFalse("derivative" in str(P_P))

  def test_issue1464(self):
    n = 6
    x = ca.SX.sym("x",n)
    u = ca.SX.sym("u")


    N = 9

    rk4 = ca.Function("f",[x,u],[x+u])

    for XX,XFunction in [(ca.SX,ca.Function),(ca.MX,ca.Function)]:

      g = []
      g2 = []


      V = XX.sym("V",(N+1)*n+N)
      VX,VU = ca.vertsplit(V,[0,(N+1)*n,(N+1)*n+N])

      VXk = ca.vertsplit(VX,n)
      VUk = ca.vertsplit(VU,1)

      for k in range(N):

          xf = rk4(VXk[k],VUk[k])

          xfp = ca.vertsplit(xf,int(n/2))
          vp = ca.vertsplit(VXk[k+1],int(n/2))

          g.append(xfp[0] - vp[0])
          g.append(xfp[1] - vp[1])

          g2.append(xf-VXk[k+1])

      for i in range(2):
        f = XFunction("nlp",[V],[ca.vertcat(*g)],{"ad_weight_sp":i})

        assert f.jac_sparsity(0, 0).nnz()==162

        f2 = XFunction("nlp",[V],[ca.vertcat(*g2)],{"ad_weight_sp":i})

        assert f2.jac_sparsity(0, 0).nnz()==162

  def test_callback(self):
    class mycallback(ca.Callback):
      def __init__(self, name, opts={}):
        ca.Callback.__init__(self)
        self.construct(name, opts)

      def eval(self,argin):
        return [argin[0]**2]

    foo = mycallback("my_f")

    x = ca.MX.sym('x')
    y = foo(x)

    f = ca.Function("f",[x],[y])

    out = f(5)

    self.checkarray(out,25)

  def test_callback_buffer(self):
    class mycallback(ca.Callback):
      def __init__(self, name, opts={}):
        ca.Callback.__init__(self)
        self.construct(name, opts)
      def has_eval_buffer(self): return True
      def eval_buffer(self, arg, res):
        a = np.frombuffer(arg[0], dtype=np.float64)
        r = np.frombuffer(res[0], dtype=np.float64)
        r[:] = a**2
        return 0

    foo = mycallback("my_f")

    x = ca.MX.sym('x')
    y = foo(x)

    f = ca.Function("f",[x],[y])

    out = f(5)

    self.checkarray(out,25)

    class mycallback(ca.Callback):
      def __init__(self, name, opts={}):
        ca.Callback.__init__(self)
        self.construct(name, opts)
      def has_eval_buffer(self): return True
      def get_n_in(self): return 3
      def get_n_out(self): return 2
      def get_sparsity_in(self, i):
        if i==0:
          return ca.Sparsity.dense(1,1)
        elif i==1:
          return ca.Sparsity.dense(3,1)
        else:
          return ca.Sparsity.dense(3,3)
      def get_sparsity_out(self, i):
        if i==0:
          return ca.Sparsity.dense(3,1)
        else:
          return ca.Sparsity.dense(3,3)
      def eval_buffer(self, arg, res):
        a = np.frombuffer(arg[0], dtype=np.float64)
        b = np.frombuffer(arg[1], dtype=np.float64)
        c = np.frombuffer(arg[2], dtype=np.float64).reshape((3,3), order='F')
        print(c)
        r0 = np.frombuffer(res[0], dtype=np.float64)
        r1 = np.frombuffer(res[1], dtype=np.float64).reshape((3,3), order='F')
        r0[:] = np.dot(a*c,b)
        r1[:,:] = c**2
        return 0

    foo = mycallback("my_f")

    a = 3
    b = ca.DM([1,2,3])
    c = ca.DM([[1,2,3],[4,5,6],[7,8,9]])

    res = foo(a,b,c)

    self.checkarray(res[0],a*c @ b)
    self.checkarray(res[1],c**2)

  def test_callback_errors(self):
    class mycallback(ca.Callback):
      def __init__(self, name, opts={}):
        ca.Callback.__init__(self)
        self.construct(name, opts)
      def eval(self,argin):
        raise Exception("foobar")

    foo = mycallback("my_f")

    x = ca.MX.sym('x')
    y = foo(x)

    f = ca.Function("f",[x],[y])

    try:
      f(3)
    except Exception as e:
      self.assertTrue("foobar" in str(e))

  def test_callback_convolution(self):
    import casadi as ca

    class MySVD(ca.Callback):
        def __init__(self, name, n, m, fd):
            ca.Callback.__init__(self)
            self.n = n
            self.m = m
            self.fd = fd
            self.k = min(n, m)
            self.fwd = None
            self.adj = None
            self.construct(name)

        def get_n_in(self):
            return 1

        def get_n_out(self):
            return 3

        def get_sparsity_in(self, i):
            return ca.Sparsity.dense(self.n, self.m)

        def get_sparsity_out(self, i):
            if i == 0:
                return ca.Sparsity.dense(self.n, self.k)
            elif i == 1:
                return ca.Sparsity.dense(self.k)
            elif i == 2:
                return ca.Sparsity.dense(self.k, self.m)
            else:
                raise ValueError('No such index')

        def eval(self, arg):
            u, s, vh = np.linalg.svd(np.array(arg[0]).reshape(self.n, self.m), 
                                      full_matrices=False)
            return [u, s, vh]
    class Convolution(ca.Callback):
        def __init__(self, A, flag_transp=False):
            ca.Callback.__init__(self)
            self.A = A
            if flag_transp:
                self.m, self.n = A.shape
            else:
                self.n, self.m = A.shape
            self.data = []
            self.flag_transp = flag_transp
            self.construct('Convolution')

        def get_transpose(self):
            out = Convolution(self.A, not self.flag_transp)
            self.data.append(out)
            return out

        def get_sparsity_in(self, i):
            return ca.Sparsity.dense(self.m, 1)

        def get_sparsity_out(self, i):
            return ca.Sparsity.dense(self.n, 1)

        def eval(self, arg):
            x = arg[0]
            # Fill in smarter numerical code
            if self.flag_transp:
                y = self.A.T @ x
            else:
                y = self.A @ x
            return [y]

        def get_n_in(self):
            return 1

        def get_n_out(self):
            return 1

        def has_forward(self, nfwd):
            return True

        def has_reverse(self, nadj):
            return True

        def get_forward(self, nfwd, name, inames, onames, opts):
            mf = self.map(nfwd, 'serial')
            fwd_in = ca.MX.sym('x', self.m, nfwd)
            fwd_out = mf(fwd_in)
            arg = ca.MX.sym('x', self.m, 1)
            dummy = ca.MX.sym('x', ca.Sparsity(self.n, 1))
            return ca.Function(name, [arg, dummy, fwd_in], [fwd_out],
                               inames, onames, opts)

        def get_reverse(self, nadj, name, inames, onames, opts):
            mf = self.get_transpose().map(nadj, 'serial')
            adj_in = ca.MX.sym('x', self.n, nadj)
            adj_out = mf(adj_in)
            arg = ca.MX.sym('x', self.m, 1)
            dummy = ca.MX(self.n, 1)
            return ca.Function(name, [arg, dummy, adj_in], [adj_out],
                               inames, onames, opts)

        def has_jacobian(self):
            return False

    n = 3
    m = 3
    foo = MySVD('foo', n, m, True)

    ca.DM.rng(1)

    x = ca.DM.rand(n,m)
    X = ca.MX.sym('x',n,m)
    Y = foo(X)

    A = ca.DM.rand(5,5)

    foo = Convolution(A, False)
    foo.forward(1)
    foo.jacobian()

  def test_mapdict(self):
    x = ca.SX.sym("x")
    y = ca.SX.sym("y",2)
    z = ca.SX.sym("z",2,2)
    v = ca.SX.sym("z",ca.Sparsity.upper(3))

    fun = ca.Function("f",{"x":x,"y":y,"z":z,"v":v,"I":z @ y+x,"II":ca.sin(y*x).T,"III":v/x},["x","y","z","v"],["I","II","III"])

    n = 2

    X = [ca.MX.sym("x") for i in range(n)]
    Y = [ca.MX.sym("y",2) for i in range(n)]
    Z = [ca.MX.sym("z",2,2) for i in range(n)]
    V = [ca.MX.sym("z",ca.Sparsity.upper(3)) for i in range(n)]

    res = fun.map(n).call({"x":ca.horzcat(*X),"y":ca.horzcat(*Y),"z":ca.horzcat(*Z),"v":ca.horzcat(*V)})

    res2 = fun.map(n).call([ca.horzcat(*X),ca.horzcat(*Y),ca.horzcat(*Z),ca.horzcat(*V)])

    F = ca.Function("F",X+Y+Z+V,res2)
    F2 = ca.Function("F",X+Y+Z+V,[res["I"],res["II"],res["III"]])

    np.random.seed(0)
    X_ = [ ca.DM(i.sparsity(),np.random.random(i.nnz())) for i in X ]
    Y_ = [ ca.DM(i.sparsity(),np.random.random(i.nnz())) for i in Y ]
    Z_ = [ ca.DM(i.sparsity(),np.random.random(i.nnz())) for i in Z ]
    V_ = [ ca.DM(i.sparsity(),np.random.random(i.nnz())) for i in V ]

    self.checkfunction(F,F2,inputs=X_+Y_+Z_+V_,jacobian=False,hessian=False,evals=False)

  @memory_heavy()
  def test_map_node(self):
    x = ca.SX.sym("x")
    y = ca.SX.sym("y",2)
    z = ca.SX.sym("z",2,2)
    v = ca.SX.sym("z",ca.Sparsity.upper(3))

    fun = ca.Function("f",[x,y,z,v],[z @ y+x,ca.sin(y*x).T,v/x])

    n = 2

    X = [ca.MX.sym("x") for i in range(n)]
    Y = [ca.MX.sym("y",2) for i in range(n)]
    Z = [ca.MX.sym("z",2,2) for i in range(n)]
    V = [ca.MX.sym("z",ca.Sparsity.upper(3)) for i in range(n)]

    for parallelization in ["serial","openmp","unroll","inline","thread"] if args.run_slow else ["serial"]:
        print(parallelization)
        res = fun.map(n, parallelization).call([ca.horzcat(*x) for x in [X,Y,Z,V]])


        F = ca.Function("F",X+Y+Z+V,list(map(ca.sin,res)))

        resref = [[] for i in range(fun.n_out())]
        for r in zip(X,Y,Z,V):
          for i,e in enumerate(map(ca.sin,fun.call(r))):
            resref[i] = resref[i] + [e]

        Fref = ca.Function("F",X+Y+Z+V,[ca.horzcat(*x) for x in resref])

        np.random.seed(0)
        X_ = [ ca.DM(i.sparsity(),np.random.random(i.nnz())) for i in X ]
        Y_ = [ ca.DM(i.sparsity(),np.random.random(i.nnz())) for i in Y ]
        Z_ = [ ca.DM(i.sparsity(),np.random.random(i.nnz())) for i in Z ]
        V_ = [ ca.DM(i.sparsity(),np.random.random(i.nnz())) for i in V ]

        for f in [F, F.expand('expand_'+F.name())]:

          self.checkfunction(f,Fref,inputs=X_+Y_+Z_+V_,sparsity_mod=args.run_slow)

  def test_map_node_light(self):
    x = ca.SX.sym("x")
    y = ca.SX.sym("y",2)
    z = ca.SX.sym("z",2,2)
    v = ca.SX.sym("z",ca.Sparsity.upper(3))

    fun = ca.Function("f",[x,y,z,v],[z @ y+x,ca.sin(y*x).T,v/x])

    n = 2

    X = [ca.MX.sym("x") for i in range(n)]
    Y = [ca.MX.sym("y",2) for i in range(n)]
    Z = [ca.MX.sym("z",2,2) for i in range(n)]
    V = [ca.MX.sym("z",ca.Sparsity.upper(3)) for i in range(n)]

    for parallelization in ["serial","openmp","unroll","inline","thread"]:
        print(parallelization)
        res = fun.map(n, parallelization).call([ca.horzcat(*x) for x in [X,Y,Z,V]])


        F = ca.Function("F",X+Y+Z+V,list(map(ca.sin,res)))

        resref = [[] for i in range(fun.n_out())]
        for r in zip(X,Y,Z,V):
          for i,e in enumerate(map(ca.sin,fun.call(r))):
            resref[i] = resref[i] + [e]

        Fref = ca.Function("F",X+Y+Z+V,[ca.horzcat(*x) for x in resref])

        np.random.seed(0)
        X_ = [ ca.DM(i.sparsity(),np.random.random(i.nnz())) for i in X ]
        Y_ = [ ca.DM(i.sparsity(),np.random.random(i.nnz())) for i in Y ]
        Z_ = [ ca.DM(i.sparsity(),np.random.random(i.nnz())) for i in Z ]
        V_ = [ ca.DM(i.sparsity(),np.random.random(i.nnz())) for i in V ]

        for f in [F, F.expand('expand_'+F.name())]:
          self.checkfunction_light(f,Fref,inputs=X_+Y_+Z_+V_,)

  def test_map_node_n_threads(self):
    x = ca.SX.sym("x")
    y = ca.SX.sym("y",2)
    z = ca.SX.sym("z",2,2)
    v = ca.SX.sym("z",ca.Sparsity.upper(3))

    fun = ca.Function("f",[x,y,z,v],[z @ y+x,ca.sin(y*x).T,v/x])

    X_ = [ ca.DM(x.sparsity(),np.random.random(x.nnz())) for i in range(10) ]
    Y_ = [ ca.DM(y.sparsity(),np.random.random(y.nnz())) for i in range(10) ]
    Z_ = [ ca.DM(z.sparsity(),np.random.random(z.nnz())) for i in range(10) ]
    V_ = [ ca.DM(v.sparsity(),np.random.random(v.nnz())) for i in range(10) ]


    print(fun.map(3,"thread",2))


    self.checkfunction_light(fun.map(2,"thread",1),fun.map(2),inputs=[ca.hcat(X_[:2]),ca.hcat(Y_[:2]),ca.hcat(Z_[:2]),ca.hcat(V_[:2])])
    self.checkfunction_light(fun.map(3,"thread",1),fun.map(3),inputs=[ca.hcat(X_[:3]),ca.hcat(Y_[:3]),ca.hcat(Z_[:3]),ca.hcat(V_[:3])])
    self.checkfunction_light(fun.map(3,"thread",2),fun.map(3),inputs=[ca.hcat(X_[:3]),ca.hcat(Y_[:3]),ca.hcat(Z_[:3]),ca.hcat(V_[:3])])
    self.checkfunction_light(fun.map(4,"thread",2),fun.map(4),inputs=[ca.hcat(X_[:4]),ca.hcat(Y_[:4]),ca.hcat(Z_[:4]),ca.hcat(V_[:4])])
    self.checkfunction_light(fun.map(4,"thread",5),fun.map(4),inputs=[ca.hcat(X_[:4]),ca.hcat(Y_[:4]),ca.hcat(Z_[:4]),ca.hcat(V_[:4])])

  @memory_heavy()
  def test_mapsum(self):
    x = ca.SX.sym("x")
    y = ca.SX.sym("y",2)
    z = ca.SX.sym("z",2,2)
    v = ca.SX.sym("z",ca.Sparsity.upper(3))

    fun = ca.Function("f",[x,y,z,v],[z @ y+x,ca.sin(y*x).T,v/x])

    n = 2

    X = [ca.MX.sym("x") for i in range(n)]
    Y = [ca.MX.sym("y",2) for i in range(n)]
    Z = [ca.MX.sym("z",2,2) for i in range(n)]
    V = [ca.MX.sym("z",ca.Sparsity.upper(3)) for i in range(n)]

    zi = 0
    for Z_alt in [Z,[ca.MX()]*3]:
      zi+= 1
      for parallelization in ["serial","openmp","unroll","thread"]:
        res = fun.mapsum([ca.horzcat(*x) for x in [X,Y,Z_alt,V]],parallelization) # Joris - clean alternative for this?

        for ad_weight_sp in [0,1]:
          F = ca.Function("F",X+Y+Z+V,list(map(ca.sin,res)),{"ad_weight": 0,"ad_weight_sp":ad_weight_sp})

          resref = [0 for i in range(fun.n_out())]  # type: list
          for r in zip(X,Y,Z_alt,V):
            for i,e in enumerate(fun.call(r)):
              resref[i] = resref[i] + e

          Fref = ca.Function("F",X+Y+Z+V,list(map(ca.sin,resref)))

          np.random.seed(0)
          X_ = [ ca.DM(i.sparsity(),np.random.random(i.nnz())) for i in X ]
          Y_ = [ ca.DM(i.sparsity(),np.random.random(i.nnz())) for i in Y ]
          Z_ = [ ca.DM(i.sparsity(),np.random.random(i.nnz())) for i in Z ]
          V_ = [ ca.DM(i.sparsity(),np.random.random(i.nnz())) for i in V ]

          inputs = X_+Y_+Z_+V_

          if parallelization!="thread":
            self.check_codegen(F,inputs=inputs)

          for f in [F,toSX_fun(F)]:
            self.checkfunction(f,Fref,inputs=inputs,sparsity_mod=args.run_slow)


  @memory_heavy()
  def test_mapsum2(self):
    x = ca.SX.sym("x")
    y = ca.SX.sym("y",2)
    z = ca.SX.sym("z",2,2)
    v = ca.SX.sym("z",ca.Sparsity.upper(3))

    fun = ca.Function("f",[x,y,z,v],[z @ y+x,ca.sin(y*x).T,v/x])

    n = 2

    X = [ca.MX.sym("x") for i in range(n)]
    Y = [ca.MX.sym("y",2) for i in range(n)]
    Z = ca.MX.sym("z",2,2)
    V = ca.MX.sym("z",ca.Sparsity.upper(3))

    for Z_alt in [Z]:

      for parallelization in ["true_map_sum","serial","openmp","unroll","thread"]:
        for ad_weight_sp in [0,1]:
          for ad_weight in [0,1]:

            if parallelization=="true_map_sum":
              F = fun.map(n,[False,False,True,True],[True,False,False],{"ad_weight_sp":ad_weight_sp,"ad_weight":ad_weight})
            else:
              F = fun.map("map",parallelization,n,[2,3],[0],{"ad_weight_sp":ad_weight_sp,"ad_weight":ad_weight})

            resref = [0 for i in range(fun.n_out())]  # type: list
            acc = 0
            bl = []
            cl = []
            for r in zip(X,Y,[Z_alt]*n,[V]*n):
              a,b,c= fun(*r)
              acc = acc + a
              bl.append(b)
              cl.append(c)

            Fref = ca.Function("F",[ca.horzcat(*X),ca.horzcat(*Y),Z,V],[acc,ca.horzcat(*bl),ca.horzcat(*cl)])

            np.random.seed(0)
            X_ = [ ca.DM(i.sparsity(),np.random.random(i.nnz())) for i in X ]
            Y_ = [ ca.DM(i.sparsity(),np.random.random(i.nnz())) for i in Y ]
            Z_ = ca.DM(Z.sparsity(),np.random.random(Z.nnz()))
            V_ = ca.DM(V.sparsity(),np.random.random(V.nnz()))

            inputs = [ca.horzcat(*X_),ca.horzcat(*Y_),Z_,V_]

            if parallelization!="thread":
              self.check_codegen(F,inputs=inputs)

            for f in [F,toSX_fun(F)]:
              self.checkfunction(f,Fref,inputs=inputs,sparsity_mod=args.run_slow)

            self.check_serialize(F,inputs=inputs)

  def test_issue4274(self):
    x=ca.MX.sym('x')

    c = ca.Function('x',[x],[x.printme(4)**2],["x"],["out"])
    local_inputs= {"x": 5}


    F = c.map(2,"thread")

    wide_inputs = {}
    for k,v in local_inputs.items():
        wide_inputs[k] = ca.horzcat(local_inputs[k],local_inputs[k]*1.2)

    print(wide_inputs)
    F(**wide_inputs)

  def test_repmatnode(self):
    x = ca.MX.sym("x",2)

    y = ca.sin(ca.repmat(x**2,1,3))

    z = ca.MX.sym("y",2,2)

    F = ca.Function("f",[x,z],[ca.sum2(ca.sum1(y))])

    x = ca.SX.sym("x",2)

    y = ca.sin(ca.repmat(x**2,1,3))
    z = ca.SX.sym("y",2,2)

    Fref = ca.Function("f",[x,z],[ca.sum2(ca.sum1(y))])

    x0 = ca.DM([1,7])
    x1 = ca.DM([[3,0],[2,4]])

    self.check_codegen(F,inputs=[x0,x1])
    self.checkfunction(F,Fref,inputs=[x0,x1])

  def test_repsumnode(self):

    x = ca.MX.sym("x",2)
    z = ca.MX.sym("y",2,2)

    F = ca.Function("f",[x,z],[ca.sin(ca.repsum((x**2).T,1,2)),(ca.cos(x**2)*2*x).T])

    x = ca.SX.sym("x",2)
    z = ca.SX.sym("y",2,2)


    Fref = ca.Function("f",[x,z],[ca.sin(ca.repsum((x**2).T,1,2)),(ca.cos(x**2)*2*x).T])

    x0 = ca.DM([1,7])
    x1 = ca.DM([[3,0],[2,4]])

    self.check_codegen(F,inputs=[x0,x1])
    self.checkfunction(F,Fref,inputs=[x0,x1])

  def test_unknown_options(self):
    x = ca.SX.sym("x")

    with self.assertRaises(Exception):
      f = SXFunction("f", [x],[x],{"fooo": False})  # pyright: ignore[reportUndefinedVariable]

    with self.assertRaises(Exception):
      f = SXFunction("f", [x],[x],{"ad_weight": "foo"})  # pyright: ignore[reportUndefinedVariable]

    if not ca.has_nlpsol("ipopt"):
      return

  @requires_nlpsol("ipopt")
  def test_unknown_options_stringvector(self):
    x = ca.SX.sym("x")
    solver = ca.nlpsol("mysolver", "ipopt", {"x":x,"f":x**2}, {"monitor": ["nlp_f"]})
    with capture_stdout() as result:
      solver = ca.nlpsol("mysolver", "ipopt", {"x":x,"f":x**2}, {"monitor": ["abc"]})
    self.assertTrue("Ignoring monitor 'abc'. Available functions: nlp_f" in result[1])

  @memory_heavy()
  def test_mapaccum(self):

    x = ca.SX.sym("x",2)
    y = ca.SX.sym("y")
    z = ca.SX.sym("z",2,2)
    v = ca.SX.sym("v",ca.Sparsity.upper(3))

    fun = ca.Function("f",[x,y,z,v],[z @ x+y,ca.sin(y*x).T,v/y])

    n = 2

    X = ca.MX.sym("x",x.sparsity())
    Y = [ca.MX.sym("y",y.sparsity()) for i in range(n)]
    Z = [ca.MX.sym("z",z.sparsity()) for i in range(n)]
    V = [ca.MX.sym("v",v.sparsity()) for i in range(n)]

    np.random.seed(0)
    X_ = ca.DM(x.sparsity(),np.random.random(x.nnz()))
    Y_ = [ ca.DM(i.sparsity(),np.random.random(i.nnz())) for i in Y ]
    Z_ = [ ca.DM(i.sparsity(),np.random.random(i.nnz())) for i in Z ]
    V_ = [ ca.DM(i.sparsity(),np.random.random(i.nnz())) for i in V ]

    for ad_weight in range(2):
      for ad_weight_sp in range(2):
        F = fun.mapaccum("map",n,[0],[0],{"ad_weight_sp":ad_weight_sp,"ad_weight": ad_weight})

        F.forward(2)

        XP = X

        Y0s = []
        Y1s = []
        Xps = []
        for k in range(n):
          XP, Y0,Y1 = fun(XP,Y[k],Z[k],V[k])
          Y0s.append(Y0)
          Y1s.append(Y1)
          Xps.append(XP)
        Fref = ca.Function("f",[X,ca.horzcat(*Y),ca.horzcat(*Z),ca.horzcat(*V)],[ca.horzcat(*Xps),ca.horzcat(*Y0s),ca.horzcat(*Y1s)])
        inputs = [X_,ca.horzcat(*Y_),ca.horzcat(*Z_),ca.horzcat(*V_)]

        for f in [F,toSX_fun(F)]:

          self.checkfunction(f,Fref,inputs=inputs)
          self.check_codegen(f,inputs=inputs)

    fun = ca.Function("f",[y,x,z,v],[z @ x+y+c.trace(v)**2,ca.sin(y*x).T,v/y])

    for ad_weight in range(2):
      for ad_weight_sp in range(2):
        F = fun.mapaccum("map",n,[1,3],[0,2],{"ad_weight_sp":ad_weight_sp,"ad_weight": ad_weight})

        XP = X
        VP = V[0]

        Y0s = []
        Y1s = []
        Xps = []
        Vps = []
        for k in range(n):
          XP, Y0, VP = fun(Y[k],XP,Z[k],VP)
          Y0s.append(Y0)
          Xps.append(XP)
          Vps.append(VP)

        Fref = ca.Function("f",[ca.horzcat(*Y),X,ca.horzcat(*Z),V[0]],[ca.horzcat(*Xps),ca.horzcat(*Y0s),ca.horzcat(*Vps)])
        inputs = [ca.horzcat(*Y_),X_,ca.horzcat(*Z_),V_[0]]

        for f in [F,toSX_fun(F)]:
          self.checkfunction(f,Fref,inputs=inputs)
          self.check_codegen(f,inputs=inputs)

  def test_mapaccum_schemes(self):

    x = ca.SX.sym("x",2)
    y = ca.SX.sym("y")
    z = ca.SX.sym("z",2,2)
    v = ca.SX.sym("v",ca.Sparsity.upper(3))

    fun = ca.Function("f",[y,z,x,v],[z @ x+y,ca.sin(y*x).T,v/y],["y","z","x","v"],["out0","out1","out2"])

    n = 2

    F = fun.mapaccum("map",n,[2],[0])

    scheme_in_fun = fun.name_in()
    scheme_out_fun = fun.name_out()

    scheme_in_F = F.name_in()
    scheme_out_F = F.name_out()

    self.assertTrue(len(scheme_in_fun),len(scheme_in_F))
    self.assertTrue(len(scheme_out_fun),len(scheme_out_F))

    for sf,sF in zip(scheme_in_fun,scheme_in_F):
      self.assertTrue(sf==sF)
    for sf,sF in zip(scheme_out_fun,scheme_out_F):
      self.assertTrue(sf==sF)

    fun = ca.Function("f",[x,y,z,v],[z @ x+y,ca.sin(y*x).T,v/y],["x","y","z","v"],["out0","out1","out2"])

    n = 2

    F = fun.mapaccum("map",n)

    self.assertTrue(len(scheme_in_fun),len(scheme_in_F))
    self.assertTrue(len(scheme_out_fun),len(scheme_out_F))

    for sf,sF in zip(scheme_in_fun,scheme_in_F):
      self.assertTrue(sf==sF)
    for sf,sF in zip(scheme_out_fun,scheme_out_F):
      self.assertTrue(sf==sF)

  # @requiresPlugin(Importer,"clang")
  # def test_jitfunction_clang(self):
  #   x = MX.sym("x")
  #   F = Function("f",[x],[x**2],{'jit':True})

  #   out = F([5])
  #   self.checkarray(out[0],25)

  # @requiresPlugin(Importer,"clang")
  # def test_clang_c(self):
  #   compiler = Importer('../data/helloworld.c', 'clang')
  #   f = external("helloworld_c", compiler)
  #   [v] = f([])
  #   self.checkarray(2.37683, v, digits=4)

  # @requiresPlugin(Importer,"clang")
  # def test_clang_cxx(self):
  #   compiler = Importer('../data/helloworld.cxx', 'clang')
  #   f = external("helloworld_cxx", compiler)
  #   [v] = f([])
  #   self.checkarray(2.37683, v, digits=4)

  # @requiresPlugin(Importer,"shell")
  # def test_shell_c(self):
  #   compiler = Importer('../data/helloworld.c', 'shell')
  #   f = external("helloworld_c", compiler)
  #   [v] = f([])
  #   self.checkarray(2.37683, v, digits=4)

  # @requiresPlugin(Importer,"shell")
  # def test_shell_cxx(self):
  #   opts = {'compiler':'g++'}
  #   compiler = Importer('../data/helloworld.cxx', 'shell', opts)
  #   f = external("helloworld_cxx", compiler)
  #   [v] = f([])
  #   self.checkarray(2.37683, v, digits=4)

  def test_depends_on(self):
    x = ca.SX.sym("x")
    y = x**2
    try:
        ca.depends_on(x,y)
    except Exception as e:
        s = str(e)
    self.assertTrue("not symbolic" in s)
    try:
        ca.Function("f",[y],[x])
    except Exception as e:
        s = str(e)
    self.assertTrue("not symbolic" in s)

  def test_1d_interpolant(self):
    grid = [[0, 1, 1.5, 2, 3]]
    values = [0, 1, 2, 5, 3]
    F = ca.interpolant('F', 'linear', grid, values)
    def same(a, b): return abs(float(a)-b)<1e-8
    pairs = [
      (3.4,3-0.4*2),
      (2.4,5-0.4*2),
      (1.6,2+3*0.1/0.5),
      (1.4,1+0.4/0.5),
      (0.4,0.4),
      (-.6,-0.6)
    ]

    X = ca.MX.sym("x")

    J = ca.Function("F",[X],[F(X)])

    for a,r in pairs:
      self.assertTrue(same(F(a), r))
      self.check_codegen(F,inputs=[a],check_serialize=True)
      self.check_serialize(F,[a])
      
    print(F.stats())

    X = ca.MX.sym("x")

    J = ca.Function("F",[X],[ca.jacobian(F(X),X)])

    pairs = [
      (3.4,-2),
      (2.4,-2),
      (1.6,6),
      (1.4,2),
      (0.4,1),
      (-.6,1),

      (1,2),
      (0.99,1),
    ]

    for a,r in pairs:
      self.assertTrue(same(J(a), r))
      self.check_codegen(J,inputs=[a])
      self.check_serialize(J,[a])

  def test_2d_interpolant(self):
    grid = [[0, 1, 4, 5],
            [0, 2, 3]]

    values = [0,   1,  8,  3,
              10, -11, 12, 13,
              20, 31, -42, 53]
    F = ca.interpolant('F', 'linear', grid, values)


    a0 = -11+0.4*(31+11)
    a1 = 12+0.4*(-42-12)
    pairs = [
      (ca.vertcat(1,2), -11),
      (ca.vertcat(1,3), 31),
      (ca.vertcat(4,2), 12),
      (ca.vertcat(4,3), -42),

      (ca.vertcat(1,2.4), a0),
      (ca.vertcat(4,2.4), a1),

      (ca.vertcat(3,2), -11+2.0/3*(12+11)),
      (ca.vertcat(3,3), 31+2.0/3*(-42-31)),

      (ca.vertcat(3,2.4), a0+2.0/3*(a1-a0))
    ]

    for a,r in pairs:
      self.checkarray(F(a), r)
      self.check_codegen(F,inputs=[a])


    X = ca.MX.sym("x",2)

    J = ca.Function("J",[X],[ca.jacobian(F(X),X)])

    jx0 = (12+11)/3.0
    jx1 = (-42-31)/3.0
    jx2 = (13-12)
    jx3 = (53+42)

    jy0 = 31+11
    jy1 = -42-12

    pairs = [
      (ca.vertcat(1,2), ca.vertcat(jx0,jy0)),
      (ca.vertcat(1,3), ca.vertcat(jx1,jy0)),
      (ca.vertcat(4,2), ca.vertcat(jx2,jy1)),
      (ca.vertcat(4,3), ca.vertcat(jx3,jy1)),

      (ca.vertcat(1,2.4), ca.vertcat(jx0+(jx1-jx0)*0.4, 31+11)),
      (ca.vertcat(4,2.4), ca.vertcat(jx2+(jx3-jx2)*0.4, -42-12)),

      (ca.vertcat(3,2), ca.vertcat(jx0,jy0+(jy1-jy0)*2.0/3)),
      (ca.vertcat(3,3), ca.vertcat(jx1,jy0+(jy1-jy0)*2.0/3)),

      (ca.vertcat(3,2.4), ca.vertcat(jx0+(jx1-jx0)*0.4,jy0+(jy1-jy0)*2.0/3)),

    ]

    for a,r in pairs:
      self.checkarray(J(a).T, r)
      self.check_codegen(J,inputs=[a])
      self.check_serialize(J,[a])

  def test_1d_interpolant_uniform(self):
    grid = [[0, 1, 2]]
    values = [0, 1, 2]
    for opts in [{"lookup_mode": ["linear"]},{"lookup_mode": ["exact"]},{"lookup_mode": ["binary"]}]:
      F = ca.interpolant('F', 'linear', grid, values, opts)
      def same(a, b): return abs(float(a)-b)<1e-8
      self.assertTrue(same(F(2.4), 2.4))
      self.assertTrue(same(F(1.4), 1.4))
      self.assertTrue(same(F(0), 0))
      self.assertTrue(same(F(1), 1))
      self.assertTrue(same(F(2), 2))
      self.assertTrue(same(F(6), 6))
      self.assertTrue(same(F(0.4), 0.4))
      self.assertTrue(same(F(-.6), -.6))

      F = ca.interpolant('F', 'linear', [np.linspace(0,1,7)], list(range(7)), {"lookup_mode": ["exact"]})

    grid = [[2, 4, 6]]
    values = [10, 7, 1]
    for opts in [{"lookup_mode": ["linear"]},{"lookup_mode": ["exact"]},{"lookup_mode": ["binary"]}]:
      F = ca.interpolant('F', 'linear', grid, values, opts)
      def same(a, b): return abs(float(a)-b)<1e-8
      self.assertTrue(same(F(1), 11.5))
      self.assertTrue(same(F(2), 10))
      self.assertTrue(same(F(3), 8.5))
      self.assertTrue(same(F(4), 7))
      self.assertTrue(same(F(5), 4))
      self.assertTrue(same(F(6), 1))
      self.assertTrue(same(F(7), -2))

      F = ca.interpolant('F', 'linear', [np.linspace(0,1,7)], list(range(7)), {"lookup_mode": ["exact"]})

  def test_2d_interpolant_uniform(self):
    grid = [[0, 1, 2], [0, 1, 2]]
    values = [0, 1, 2, 10, 11, 12, 20, 21, 22]
    for opts in [{"lookup_mode": ["linear","linear"]},{"lookup_mode": ["exact","exact"]},{"lookup_mode": ["binary","binary"]}]:
      F = ca.interpolant('F', 'linear', grid, values, opts)
      def same(a, b): return abs(float(a)-b)<1e-8
      self.assertTrue(same(F([2.4, 0.5]), 7.4))
      self.assertTrue(same(F([1.4, 0.5]), 6.4))
      self.assertTrue(same(F([0.4, 0.5]), 5.4))
      self.assertTrue(same(F([1, 0.5]), 6))
      self.assertTrue(same(F([1, 0]), 1))
      self.assertTrue(same(F([0, 0]), 0))
      self.assertTrue(same(F([0.4, 1]), 10.4))
      self.assertTrue(same(F([-.6, 0.5]), 4.4))
      self.assertTrue(same(F([-.6, 1.5]), 14.4))
      self.assertTrue(same(F([-.6, 2.5]), 24.4))
      self.assertTrue(same(F([-.6, 3.5]), 34.4))

  @skip(not scipy_interpolate)
  def test_nd_linear(self):

    N = 4
    for n_dim in range(1,5):
      grid = []
      x0 = []
      for i in range(n_dim):
        g = sorted(list(np.random.random(N)))
        grid.append(g)
        x0.append(np.mean(g))

      D = np.random.random([N]*n_dim)

      xref = scipy.interpolate.interpn(grid, D, x0)

      d = D.ravel(order='F')

      F = ca.interpolant('F', 'linear', grid, d)

      self.checkarray(F(x0),xref)


  @skip(not scipy_interpolate)
  def test_2d_bspline(self):
    import scipy.interpolate
    np.random.seed(0)

    d_knots = [list(np.linspace(0,1,5)),list(np.linspace(0,1,6))]

    data = np.random.random([len(e) for e in d_knots])
    r = np.meshgrid(*d_knots,indexing='ij')

    xyz = np.vstack(list(e.ravel(order='F') for e in r)).ravel(order='F')

    d_flat = data.ravel(order='F')

    LUT = ca.interpolant('name','bspline',d_knots,d_flat)
    LUTJ = jacobian_old(LUT, 0, 0)
    LUTH = hessian_old(LUT, 0, 0)

    self.check_codegen(LUT, [ca.vertcat(0.2,0.3)])
    self.check_serialize(LUT, [ca.vertcat(0.2,0.3)])
    #scipy.interpolate.interpn(d_knots, data, [0.2,0.3], method='splinef2d')

    interp = scipy.interpolate.RectBivariateSpline(d_knots[0], d_knots[1], data)
    for x in [0,0.01,0.1,0.2,0.9,0.99,1]:
      for y in [0,0.01,0.1,0.2,0.9,0.99,1]:
        m = LUT([x,y])
        r = interp.ev(x,y)
        self.checkarray(m,r)

        m = LUTJ([x,y])[0]
        try:
          r = [interp.ev(x,y, 1, 0), interp.ev(x,y, 0, 1)]
        except:
          r = None
        if r is not None:
          self.checkarray(m,r)

        m = LUTH([x,y])[0]
        try:
          r = ca.blockcat([[interp.ev(x,y, 2, 0),interp.ev(x,y, 1, 1)],[interp.ev(x,y, 1, 1), interp.ev(x,y, 0, 2)]])
        except:
          r = None
        if r is not None:
          self.checkarray(m,r)

  def test_1d_bspline_call_fun(self):
  
    fs = []
    for X in [ca.SX,ca.MX]:

        x = X.sym("x")
        np.random.seed(0)

        d_knots = [list(np.linspace(0,1,5))]

        data = np.random.random([len(e) for e in d_knots])
        r = np.array(d_knots)

        xyz = np.vstack(list(e.ravel(order='F') for e in r)).ravel(order='F')

        d_flat = data.ravel(order='F')

        LUT = ca.interpolant('name','bspline',d_knots,d_flat)

        f = ca.Function('f',[x],[2*LUT(ca.sin(x))])
        fs.append(f)
    self.checkfunction(fs[0],fs[1],inputs=[0.32])
        

  @skip(not scipy_interpolate)
  def test_1d_bspline(self):
    import scipy.interpolate
    np.random.seed(0)

    d_knots = [list(np.linspace(0,1,5))]

    data = np.random.random([len(e) for e in d_knots])
    r = np.array(d_knots)

    xyz = np.vstack(list(e.ravel(order='F') for e in r)).ravel(order='F')

    d_flat = data.ravel(order='F')

    LUT = ca.interpolant('name','bspline',d_knots,d_flat)
    self.check_codegen(LUT, [0.2])
    LUTJ = jacobian_old(LUT, 0, 0)
    LUTH = hessian_old(LUT, 0, 0)

    interp = scipy.interpolate.InterpolatedUnivariateSpline(d_knots[0], data)
    for x in [0,0.01,0.1,0.2,0.9,0.99,1]:
      m = LUT(x)
      r = interp(x)
      self.checkarray(m,r)

      m = LUTJ(x)[0]
      try:
        r = interp(x, 1)
      except:
        r = None
      if r is not None:
        self.checkarray(m,r)

      m = LUTH(x)[0]
      try:
        r = interp(x, 2)
      except:
        r = None
      if r is not None:
        self.checkarray(m,r)

  def test_Callback_Jacobian(self):
    x = ca.MX.sym("x")
    y = ca.MX.sym("y")

    num_inputs = [0.2,0.7]

    g = ca.Function("g", [x,y],[ca.sin(x+3*y)])

    class Fun(ca.Callback):
        # sin(x+3*y)

        def __init__(self):
          ca.Callback.__init__(self)
          self.construct("Fun", {})
        def get_n_in(self): return 2
        def get_n_out(self): return 1

        def eval(self,arg):
          x = arg[0]
          y = arg[1]
          z0 = 3*y
          z1 = x+z0
          z2 = ca.sin(z1)
          return [z2]

        def has_forward(self,nfwd): return False
        def has_reverse(self,nadj): return False

        def has_jacobian(self): return True

        def get_jacobian(self, name, inames, onames, opts):
          x = ca.SX.sym("x")
          y = ca.SX.sym("y")
          out_g = ca.SX.sym('out_g', ca.Sparsity(1,1))
          J = ca.Function(name, [x,y, out_g], [ca.cos(x+3*y), 3*ca.cos(x+3*y)], inames, onames, opts)
          return J

    f = Fun()

    self.checkfunction(f,g,inputs=num_inputs,fwd=False,adj=False,indirect=False)


  def test_Callback_errors(self):

    class Fun(ca.Callback):

        def __init__(self):
          ca.Callback.__init__(self)
          self.construct("Fun", {})
        def get_n_in(self): return 2
        def get_n_out(self): return 1

        def get_sparsity_in(i):  # pyright: ignore[reportIncompatibleMethodOverride,reportSelfClsParameterName]
          # Deliberately missing `self` and returning int -- this test
          # asserts Callback reports the error at construction time.
          return 4

        def eval(self,arg):
          x = arg[0]
          y = arg[1]

          z0 = 3*y
          z1 = x+z0
          z2 = ca.sin(z1)
          return [z2]

    try:
      f = Fun()
    except Exception as e:
      s = str(e)
      print(s)
    self.assertTrue("get_sparsity_in" in s)

  def test_Callback(self):

    x = ca.MX.sym("x")
    y = ca.MX.sym("y")

    num_inputs = [0.2,0.7]

    g = ca.Function("g", [x,y],[ca.sin(x+3*y)])

    # Simple syntax
    def getP(indirect=True):

      class Fun(ca.Callback):

        def __init__(self):
          ca.Callback.__init__(self)
          self.construct("Fun", {})

        def get_n_in(self): return 2
        def get_n_out(self): return 1

        def eval(self,arg):
          x = arg[0]
          y = arg[1]

          z0 = 3*y
          z1 = x+z0
          z2 = ca.sin(z1)
          return [z2]

      self.cb = Fun()  # pyright: ignore[reportAttributeAccessIssue]

      if not indirect:
        return self.cb  # pyright: ignore[reportAttributeAccessIssue]

      f = ca.Function("f", [x,y],[self.cb(x,y)])  # pyright: ignore[reportAttributeAccessIssue]

      return f

    for indirect in [True,False]:
      f = getP(indirect=indirect)
      self.checkfunction(f,g,inputs=num_inputs,sens_der=False,jacobian=False,gradient=False,hessian=False,evals=False)

      with self.assertRaises(Exception):
        f.gradient()  # pyright: ignore[reportAttributeAccessIssue]

      with self.assertRaises(Exception):
        jacobian_old(f, 0, 0)

      with self.assertRaises(Exception):
        f.forward(1)
      with self.assertRaises(Exception):
        f.reverse(1)

  def test_Callback_dimcheck(self):
      class Fun(ca.Callback):
        def __init__(self):
          ca.Callback.__init__(self)
          self.construct("Fun")
        def get_n_in(self): return 2
        def get_n_out(self): return 1

        def eval(self,arg):
          return [2, 1]
      f = Fun()

      s = ""
      try:
        f(2)
      except Exception as e:
        s = str(e)
      self.assertTrue("Incorrect number of inputs" in s)
      class Fun(ca.Callback):
        def __init__(self):
          ca.Callback.__init__(self)
          self.construct("Fun")
        def get_n_in(self): return 2
        def get_n_out(self): return 1

        def eval(self,arg):
          return [2, 1]
      f = Fun()

      s = ""
      try:
        f(2,3)
      except Exception as e:
        s = str(e)
      self.assertTrue("Expected 1 output" in s)
      s = ""
      class Fun(ca.Callback):
        def __init__(self):
          ca.Callback.__init__(self)
          self.construct("Fun")
        def get_n_in(self): return 2
        def get_n_out(self): return 1

        def eval(self,arg):
          return [ca.DM.zeros(2,2)]
      f = Fun()
      try:
        f(2,3)
      except Exception as e:
        s = str(e)
      self.assertTrue("Shape mismatch" in s)

  def test_Callback_sens(self):
    x = ca.MX.sym("x")
    y = ca.MX.sym("y")

    num_inputs = [0.2,0.7]

    g = ca.Function("g", [x,y],[ca.sin(x+3*y)])

    def getP(has_fwd=True,has_adj=True,indirect=True):

      class Fun(ca.Callback):
        # sin(x+3*y)

        def __init__(self,opts):
          ca.Callback.__init__(self)
          self.construct("Fun", opts)
        def get_n_in(self): return 2
        def get_n_out(self): return 1

        def eval(self,arg):
          x = arg[0]
          y = arg[1]
          z0 = 3*y
          z1 = x+z0
          z2 = ca.sin(z1)
          return [z2]

        if has_fwd:
          def has_forward(self,nfwd): return nfwd==1
          def get_forward(self,nfwd,name,inames,onames,opts):
            assert(nfwd==1)
            class ForwardFun(ca.Callback):
              # sin(x+3*y)

              def __init__(self):
                ca.Callback.__init__(self)
                self.construct(name, {"verbose":True})
              def get_n_in(self): return 2+1+2
              def get_n_out(self): return 1

              def eval(self,arg):
                x,y = arg[0],arg[1]
                z = arg[2]
                seeds = arg[3:]

                z0 = 3*y
                z1 = x+z0
                z2 = ca.sin(z1)

                ret = []

                for i in range(3,len(arg),2):
                  dx = arg[i]
                  dy = arg[i+1]
                  dz0 = 3*dy
                  dz1 = dx+dz0
                  dz2 = ca.cos(z1)*dz1
                  ret.append(dz2)

                return ret

            self.cb_fwd = ForwardFun()
            return self.cb_fwd

        if has_adj:
          def has_reverse(self,nadj): return nadj==1
          def get_reverse(self,nadj,name,inames,onames,opts):
            assert(nadj==1)
            class BackwardFun(ca.Callback):
              # sin(x+3*y)

              def __init__(self):
                ca.Callback.__init__(self)
                self.construct(name, {"verbose":True})
              def get_n_in(self): return 2+1+1
              def get_n_out(self): return 2

              def eval(self,arg):
                x,y = arg[0],arg[1]
                z = arg[2]
                seeds = arg[3:]

                z0 = 3*y
                z1 = x+z0
                z2 = ca.sin(z1)

                ret = []

                for i in range(3,len(arg)):
                  z_bar = arg[i]
                  bx = 0
                  by = 0
                  bz1 = 0
                  bz0 = 0

                  bz2 = z_bar
                  bz1 += bz2*ca.cos(z1)
                  bx+= bz1;bz0+= bz1
                  by+= 3*bz0
                  ret.append(bx)
                  ret.append(by)
                return ret

            self.cb_rev = BackwardFun()
            return self.cb_rev

      opts = {"verbose":True}
      self.cb = Fun(opts)  # pyright: ignore[reportAttributeAccessIssue]
      f = self.cb  # pyright: ignore[reportAttributeAccessIssue]

      if not indirect:
        return f

      f = ca.Function("f", [x,y],[f(x,y)])

      return f

    for indirect in [True,False]:
      f = getP(has_fwd=True,has_adj=True,indirect=indirect)

      self.checkfunction(f,g,inputs=num_inputs,sens_der=False,hessian=False,evals=1)

      f = getP(has_fwd=True,has_adj=False,indirect=indirect)

      self.checkfunction(f,g,inputs=num_inputs,sens_der=False,hessian=False,adj=False,evals=1)

      f = getP(has_fwd=False,has_adj=True,indirect=indirect)

      self.checkfunction(f,g,inputs=num_inputs,sens_der=False,hessian=False,fwd=False,evals=1)

      f = getP(has_fwd=True,has_adj=False,indirect=indirect)

      self.checkfunction(f,g,inputs=num_inputs,sens_der=False,hessian=False,adj=False,evals=1)

      f = getP(has_fwd=False,has_adj=True,indirect=indirect)

      self.checkfunction(f,g,inputs=num_inputs,sens_der=False,hessian=False,fwd=False,evals=1)

  def test_callback_jacobian_sparsity(self):

    x = ca.MX.sym("x",2)

    h = 1e-7
    for with_jac_sparsity in [True, False]:

      calls = []

      class Fun(ca.Callback):

          def __init__(self):
            ca.Callback.__init__(self)
            self.construct("Fun", {"enable_fd":True,"fd_method":"forward","fd_options":{"h": h,"h_iter":False}})
          def get_n_in(self): return 1
          def get_n_out(self): return 1
          def get_sparsity_in(self,i): return ca.Sparsity.dense(2,1)
          def get_sparsity_out(self,i): return ca.Sparsity.dense(2,1)

          def eval(self,arg):
            x = arg[0]
            calls.append((x-ca.DM([5,7]))/h)
            return [x**2]

          def has_forward(self,nfwd): return False
          def has_reverse(self,nadj): return False

          def has_jac_sparsity(self,oind,iind): return with_jac_sparsity

          def get_jac_sparsity(self,oind,iind,symmetric):
            assert(oind == 0)
            assert(iind == 0)
            return ca.Sparsity.diag(2)

      f = Fun()
      y = ca.jacobian(f(x),x)

      F = ca.Function('F',[x],[y])

      with capture_stdout() as out:
       J = F([5,7])
      calls = ca.hcat(calls)
      J_ref = ca.DM([[5*2,0],[0,7*2]])
      if with_jac_sparsity:
        self.checkarray(calls,ca.DM([[0,0],[1,1]]).T,digits=5)
        J_ref = ca.sparsify(J_ref)
      else:
        self.checkarray(calls,ca.DM([[0,0],[1,0],[0,1]]).T,digits=5)

      self.checkarray(J,J_ref,digits=5)
      self.assertTrue(J.sparsity()==J_ref.sparsity())


  @requires_nlpsol("ipopt")
  def test_common_specific_options(self):

      x = ca.SX.sym("x")

      nlp = {"x": x, "f": x**4}

      with capture_stdout() as out:
        solver = ca.nlpsol("solver","ipopt",nlp)
      self.assertTrue("nlp_f" not in out[0])
      solver = ca.nlpsol("solver","ipopt",nlp,{"common_options":{"final_options" : {"print_in":True}}})
      with capture_stdout() as out:
        solver(x0=0.1)
      print(out[0])
      self.assertTrue("Function nlp_f" in out[0])
      self.assertTrue("Function nlp_grad_f" in out[0])
      solver = ca.nlpsol("solver","ipopt",nlp,{"specific_options":{ "nlp_f" : {"final_options" : {"print_in":True}}}})
      with capture_stdout() as out:
        solver(x0=0.1)
      self.assertTrue("Function nlp_f" in out[0])
      self.assertFalse("Function nlp_grad_f" in out[0])
      solver = ca.nlpsol("solver","ipopt",nlp,{"common_options":{"final_options" : {"print_in":True}},"specific_options":{ "nlp_f" : {"final_options" : {"print_in":False}}}})
      with capture_stdout() as out:
        solver(x0=0.1)
      self.assertFalse("Function nlp_f" in out[0])
      self.assertTrue("Function nlp_grad_f" in out[0])
      solver = ca.nlpsol("solver","ipopt",nlp,{"common_options":{"final_options" : {"print_in":False}},"specific_options":{ "nlp_f" : {"final_options" : {"print_in":True}}}})
      with capture_stdout() as out:
        solver(x0=0.1)
      self.assertTrue("Function nlp_f" in out[0])
      self.assertFalse("Function nlp_grad_f" in out[0])

      with capture_stdout() as out:
        solver = ca.nlpsol("solver","ipopt",nlp)
      self.assertTrue(len(out[1])==0)
      with capture_stdout() as out:
        solver = ca.nlpsol("solver","ipopt",nlp,{"specific_options":{ "nlp_foo" : {"helper_options" : {"verbose":True}}}})
      self.assertTrue("Ignoring" + out[1])
      self.assertTrue("nlp_g" in out[1])
      with self.assertRaises(Exception):
        solver = ca.nlpsol("solver","ipopt",nlp,{"specific_options":{ "nlp_foo" : 3}})


  @requires_nlpsol("ipopt")
  def test_oracle_options(self):
    ca.DM.set_precision(16)

    N = 5
    x = ca.MX.sym("x",N,1)

    options = {}
    options["ipopt"] = {"hessian_approximation":"limited-memory"}
    options["oracle_options"] = {"enable_fd":True,"enable_forward":False,"enable_reverse":False,"print_in":True,"fd_method":"central"}

    solver = ca.nlpsol("solver","ipopt",{"x":x,"f":ca.dot(x,x)},options)

    with capture_stdout() as out:
      solver(x0=1)


    ca.DM.set_precision(6)

    self.assertTrue("[1, 1, 1, 1, 1]" in out[0])


  @requires_expm("slicot")
  @memory_heavy()
  def test_expm(self):
      eps = 1e-6
      t = ca.MX.sym('t')
      tnum = 0.2

      n = 4

      np.random.seed(0)
      Anum = np.random.random((n,n))
      Bnum = np.random.random((n,2))
      Bb = np.random.random((n,2))

      dA = np.random.random((n,n))
      Yb = np.random.random((n,2))

      def expm(A):

        n = A.shape[0]
        x = ca.MX.sym('x',n)
        As = ca.MX.sym('A',n,n)

        dae = {'x':x,'p':ca.vec(As),'ode':As @ x}
        intg = ca.integrator('intg','cvodes',dae,{'reltol':1e-14,'abstol':1e-14})

        Intg = intg.map('identity','serial',n,[1],[])

        out = Intg(x0=ca.DM.eye(n),p=ca.vec(As))
        expmF = ca.Function('expm',[As],[out["xf"]])
        return expmF(A)

      A = ca.MX.sym("A",n,n)
      t = ca.MX.sym("t")
      fr = ca.Function('fr',[A,t],[expm(A*t)])
      f = ca.Function('f',[A,t],[ca.expm(A*t)])

      self.checkfunction(fr,f,inputs=[Anum, 1.1],digits=8)

      fr = ca.Function('fr',[t],[expm(Anum*t)])
      f = ca.Function('f',[t],[ca.expm_const(Anum,t)])

      self.checkfunction(fr,f,inputs=[1.1],digits=8)

      JA = ca.jacobian(ca.expm(A*t),A)
      Jt = ca.jacobian(ca.expm(A*t),t)

      self.assertTrue(JA.nnz()==n**4)
      self.assertTrue(Jt.nnz()==n**2)

      JA = ca.jacobian(ca.expm_const(A,t),A)
      Jt = ca.jacobian(ca.expm_const(A,t),t)

      self.assertTrue(JA.nnz()==0)
      self.assertTrue(Jt.nnz()==n**2)

  @memory_heavy()
  def test_conditional(self):

    np.random.seed(5)

    x = ca.MX.sym('x',2,2)
    y = ca.MX.sym('y',2,2)

    sp1 = ca.MX.sym('y',ca.Sparsity.lower(2))
    sp2 = ca.MX.sym('z',ca.Sparsity.diag(2))

    f1 = ca.Function("f",[sp2,x],[x**2,x*sp2])
    f2 = ca.Function("f",[sp1,x],[2*x**2,ca.sin(sp1)])
    f3 = ca.Function("f",[sp1,sp2],[sp1*sp2,sp1+sp2])

    F = ca.Function.conditional("test",[f1,f2], f3)
    Fsx = F.expand()

    A = np.random.random((2,2))
    B = np.random.random((2,2))

    for i in range(-1,3):
      self.checkfunction(F,Fsx,inputs = [i,A,B])
      self.check_codegen(F,inputs=[i,A,B])

  def test_max_num_dir(self):
    x = ca.MX.sym("x",10)

    f = ca.Function("ffff",[x],[ca.DM.ones(10,10) @ x],{"max_num_dir":4,"verbose":True})
    f = f.expand()


    y = f.call([ca.sin(x)],False,True)[0]



    acc = ca.Function('acc',[x],[y])
    acc = acc.mapaccum('accd',5)


    cons = ca.vec(acc(x))

    os.system("callgrind_control -i on")

    e = ca.jacobian(cons,x)
    os.system("callgrind_control -i off")

    g = ca.Function('g',[x],[e])


    c = ca.CodeGenerator('me')
    c.add(g)
    code= c.dump()

    self.assertTrue("fwd4_ffff" in code)

  def test_max_num_dir(self):
    x = ca.MX.sym("x")

    f = ca.Function("ffff",[x],[ca.sin(x)])
    fm = f.mapaccum("mapaccum",100,1,{"base":4})

    c = ca.CodeGenerator('me')
    c.add(fm)
    code= c.dump()

    self.assertTrue("ffff_acc4_acc4_acc4" in code)
    
    
  def test_codegen_with_jac_sparsity(self):
  
    if not args.run_slow: return
    x = ca.MX.sym("x",3)
    y = ca.MX.sym("y",3,3)
    f = ca.Function("f",[x,y],[x**3,y @ x])
    c = ca.CodeGenerator('me')
    c.add(f, True)
    
    fJ = f.jacobian()
    
    ca.DM.rng(0)
    x0 = ca.DM.rand(x.sparsity())
    y0 = ca.DM.rand(y.sparsity())
    cg = self.check_codegen(f,inputs=[x0,y0],with_jac_sparsity=True,external_opts={"enable_fd":True})
    FJ = cg["F"].jacobian()
    for i in range(fJ.n_out()):
        self.assertEqual(fJ.sparsity_out(i),FJ.sparsity_out(i))
    cg = self.check_codegen(f,inputs=[x0,y0],with_jac_sparsity=False,external_opts={"enable_fd":True})
    FJ = cg["F"].jacobian()
    for i in range(fJ.n_out()):
        self.assertTrue(FJ.sparsity_out(i).is_dense())

  def test_2d_linear_multiout(self):
    np.random.seed(0)

    d_knots = [list(np.linspace(0,1,5)),list(np.linspace(0,1,6))]

    data0 = np.random.random([len(e) for e in d_knots])
    data1 = np.random.random([len(e) for e in d_knots])
    r = np.meshgrid(*d_knots,indexing='ij')

    xyz = np.vstack(list(e.ravel(order='F') for e in r)).ravel(order='F')

    d_flat0 = data0.ravel(order='F')

    LUT0 = ca.interpolant('name','linear',d_knots,d_flat0)

    d_flat1 = data1.ravel(order='F')

    LUT1 = ca.interpolant('name','linear',d_knots,d_flat1)

    data = np.vstack((data0.ravel(order='F'),data1.ravel(order='F'))).ravel(order='F')

    d_flat = data.ravel(order='F')


    LUT = ca.interpolant('name','linear',d_knots,d_flat)


    x = ca.MX.sym("x")
    y = ca.MX.sym("y")


    LUT_sep = ca.Function('f',[x,y],[ca.vertcat(LUT0(ca.vertcat(x,y)),LUT1(ca.vertcat(x,y)))])
    LUT = ca.Function('f',[x,y],[LUT(ca.vertcat(x,y))])

    self.checkfunction(LUT,LUT_sep, inputs=[0.2,0.333])
    self.check_codegen(LUT,inputs=[0.2,0.333])
    self.check_serialize(LUT,inputs=[0.2,0.333])

    LUT_param = ca.interpolant('name','linear',d_knots,2)
    f = ca.Function('LUTp',[x,y],[LUT_param(ca.vertcat(x,y),d_flat)])
    self.checkfunction(LUT,f, inputs=[0.2,0.333])
    self.check_codegen(f,inputs=[0.2,0.333])
    self.check_serialize(f,inputs=[0.2,0.333])

    d_knots_cat = ca.vcat(d_knots[0]+d_knots[1])

    LUT_param = ca.interpolant('name','linear',[5,6],2)
    f = ca.Function('LUTp',[x,y],[LUT_param(ca.vertcat(x,y),d_knots_cat,d_flat)])
    self.checkfunction(LUT,f, inputs=[0.2,0.333])
    self.check_codegen(f,inputs=[0.2,0.333])
    self.check_serialize(f,inputs=[0.2,0.333])

    LUT_param = ca.interpolant('name','linear',[5,6],d_flat)
    f = ca.Function('LUTp',[x,y],[LUT_param(ca.vertcat(x,y),d_knots_cat)])
    self.checkfunction(LUT,f, inputs=[0.2,0.333])
    self.check_codegen(f,inputs=[0.2,0.333])
    self.check_serialize(f,inputs=[0.2,0.333])

  def test_2d_bspline_multiout(self):
    np.random.seed(0)

    d_knots = [list(np.linspace(0,1,5)),list(np.linspace(0,1,6))]

    data0 = np.random.random([len(e) for e in d_knots])
    data1 = np.random.random([len(e) for e in d_knots])
    r = np.meshgrid(*d_knots,indexing='ij')

    xyz = np.vstack(list(e.ravel(order='F') for e in r)).ravel(order='F')

    d_flat0 = data0.ravel(order='F')

    LUT0 = ca.interpolant('name','bspline',d_knots,d_flat0)

    d_flat1 = data1.ravel(order='F')

    LUT1 = ca.interpolant('name','bspline',d_knots,d_flat1)

    data = np.vstack((data0.ravel(order='F'),data1.ravel(order='F'))).ravel(order='F')

    d_flat = data.ravel(order='F')


    LUT = ca.interpolant('name','bspline',d_knots,d_flat)


    x = ca.MX.sym("x")
    y = ca.MX.sym("y")


    LUT_sep = ca.Function('f',[x,y],[ca.vertcat(LUT0(ca.vertcat(x,y)),LUT1(ca.vertcat(x,y)))])
    LUT = ca.Function('f',[x,y],[LUT(ca.vertcat(x,y))])

    self.checkfunction(LUT,LUT_sep, inputs=[0.2,0.333])
    self.check_codegen(LUT,inputs=[0.2,0.333])
    self.check_serialize(LUT,inputs=[0.2,0.333])

    LUT_param = ca.interpolant('name','bspline',d_knots,2)
    f = ca.Function('LUTp',[x,y],[LUT_param(ca.vertcat(x,y),d_flat)])
    self.checkfunction(LUT,f, inputs=[0.2,0.333])
    self.check_codegen(f,inputs=[0.2,0.333])
    self.check_serialize(f,inputs=[0.2,0.333])


  def test_parametric_bspline(self):
    knots = [[0,0,0,0,0.2,0.5,0.8,1,1,1,1],[0,0,0,0.1,0.5,0.9,1,1,1]]
    x=ca.MX.sym("x",2)
    data = np.random.random((7,6,2)).ravel(order='F')
    y = ca.bspline(x,data,knots,[3,2],2)

    C = ca.MX.sym("C",data.shape[0],1)
    Y = ca.bspline(x,C,knots,[3,2],2)
    f = ca.Function('f',[x],[y])
    F = ca.Function('f',[x,C],[Y])
    F = ca.Function('f',[x],[F(x,data)])
    self.checkfunction(f,F,inputs=[ca.vertcat(0.3,0.4)])
    self.check_codegen(F,inputs=[ca.vertcat(0.3,0.4)])
    self.check_serialize(F,inputs=[ca.vertcat(0.3,0.4)])

    # Parametric knots (MX knots, always inlined)
    K = [ca.MX.sym("K0",len(knots[0])), ca.MX.sym("K1",len(knots[1]))]
    Y_pk = ca.bspline(x,C,K,[3,2],2)
    # Substitute fixed knots to compare with reference
    knots_dm = [ca.DM(k) for k in knots]
    Y_sub = ca.substitute([Y_pk], K, knots_dm)[0]
    F_pk = ca.Function('f_pk',[x,C],[Y_sub])
    F_pk = ca.Function('f_pk',[x],[F_pk(x,data)])
    self.checkfunction(f,F_pk,inputs=[ca.vertcat(0.3,0.4)])
    self.check_codegen(F_pk,inputs=[ca.vertcat(0.3,0.4)],std="c99")


  def test_smooth_linear(self):
    np.random.seed(0)

    d_knots = [list(np.linspace(0,1,5)),list(np.linspace(0,1,6))]

    data0 = np.random.random([len(e) for e in d_knots])
    data1 = np.random.random([len(e) for e in d_knots])
    r = np.meshgrid(*d_knots,indexing='ij')

    data = np.vstack((data0.ravel(order='F'),data1.ravel(order='F'))).ravel(order='F')

    d_flat = data.ravel(order='F')


    LUT_linear = ca.interpolant('name','linear',d_knots,d_flat)


    LUT = ca.interpolant('name','bspline',d_knots,d_flat,{"algorithm": "smooth_linear","smooth_linear_frac":0.1})
    LUT2 = ca.interpolant('name','bspline',d_knots,d_flat,{"algorithm": "smooth_linear","smooth_linear_frac":0.05})


    # Far away from points: almost identical
    diff1 = float(ca.norm_1(LUT([0.2,0.333])-LUT_linear([0.2,0.333])))
    # Large difference near the knots
    diff2 = float(ca.norm_2(LUT([0.25,0.4])-LUT_linear([0.25,0.4])))
    diff22 = float(ca.norm_2(LUT2([0.25,0.4])-LUT_linear([0.25,0.4])))

    self.assertTrue(diff1<=1e-10)
    self.assertTrue(diff2>=100*diff1)
    self.assertTrue(diff22<=diff2*0.51)
    self.assertTrue(diff22>=diff2*0.49)

    self.checkarray(LUT2([0.26,0.39]),ca.DM([0.9261362392504342, 0.9157189108791507]))


    xy = ca.MX.sym("xy",2)
    LUT_param = ca.interpolant('name','bspline',d_knots,2,{"algorithm": "smooth_linear","smooth_linear_frac":0.1})
    f = ca.Function('LUTp',[xy],[LUT_param(xy,d_flat)])
    self.checkfunction(LUT,f, inputs=[ca.vertcat(0.2,0.333)])
    self.check_codegen(f,inputs=[ca.vertcat(0.2,0.333)],main=True)
    self.check_serialize(f,inputs=[ca.vertcat(0.2,0.333)])


  def test_codegen_avoid_stack(self):
    x = ca.SX.sym("x",3,3)
    f = ca.Function('f',[x],[ca.det(x)])
    np.random.seed(0)
    self.check_codegen(f,inputs=[np.random.random((3,3))])
    self.check_codegen(f,inputs=[np.random.random((3,3))], opts={"avoid_stack": True})


  @memory_heavy()
  def test_serialize(self):
    for opts in [{"debug":True},{}]:
      x = ca.SX.sym("x")
      y = x+3
      z = ca.sin(y)

      f = ca.Function('f',[x],[z])
      fs = ca.Function.deserialize(f.serialize(opts))

      self.checkfunction(f,fs,inputs=[2])

      x = ca.SX.sym("x")
      y = x+3
      z = ca.sin(y)

      f = ca.Function('f',[x],[z,np.nan,-np.inf,np.inf])
      fs = ca.Function.deserialize(f.serialize(opts))
      self.checkfunction(f,fs,inputs=[2])

      x = ca.SX.sym("x")
      y = ca.SX.sym("y", ca.Sparsity.lower(3))
      z = x+y
      z1 = ca.sparsify(ca.vertcat(z[0],0,z[1]))
      z2 = z.T

      f = ca.Function('f',[x,y],[z1,z2,x**2],["x","y"],["a","b","c"])
      fs = ca.Function.deserialize(f.serialize(opts))

      self.assertEqual(fs.name_in(0), "x")
      self.assertEqual(fs.name_out(0), "a")
      self.assertEqual(fs.name(), "f")

      self.checkfunction(f,fs,inputs=[3.7,np.array([[1,0,0],[2,3,0],[4,5,6]])],hessian=False)


      fs = ca.Function.deserialize(f.serialize(opts))
      self.checkfunction(f,fs,inputs=[3.7,np.array([[1,0,0],[2,3,0],[4,5,6]])],hessian=False)

      x = ca.SX.sym("x")
      p = ca.SX.sym("p")

      f = ca.Function('f',[x],[p], {"allow_free":True})

      #SXFunction with free parameters
      pickle.loads(pickle.dumps(f))


      x = ca.MX.sym("x")
      f = ca.Function('f',[x],[x**2])

      fs = ca.Function.deserialize(f.serialize(opts))
      self.checkfunction(f,fs,inputs=[3.7],hessian=False)


      x = ca.MX.sym("x")
      y = ca.MX.sym("y",2)

      w = ca.if_else(x, ca.atan2(3*ca.norm_fro(y)*y,x), x-y, True)
      z = ca.sin(2*x)*w[0]+1
      g = ca.Function("g",[x,y],[w-x])
      gmap = g.map(2, "thread", 2)
      gmapsx = gmap.expand()

      q = gmap(ca.horzcat(2*x,x-y[1]),ca.horzcat(z+y,ca.cos(z+y)))+1/gmapsx(ca.horzcat(2*x,x-y[1]),ca.repmat(z+y,1,2))

      q = ca.solve(q,2*y, "qr")
      q+= ca.bilin(ca.DM([[1,3],[7,8]]),q,2*q)

      f = ca.Function("f",[x,y],[q+1,ca.jacobian(q, ca.vertcat(x, y))])

      fs = ca.Function.deserialize(f.serialize(opts))
      self.checkfunction(f,fs,inputs=[1.1, ca.vertcat(2.7,3)],hessian=False)

  @memory_heavy()
  def test_serialize_recursion_limit(self):
      for X in [ca.SX,ca.MX]:
        x = X.sym("x")

        y = x
        for i in range(10000):
          y = ca.sin(y)

        f = ca.Function('foo',[x],[y])
        ca.Function.deserialize(f.serialize())


  def test_string(self):
    x=ca.MX.sym("x")
    y=ca.MX.sym("y")
    f = ca.Function('f',[x],[],["x"],[])
    self.assertTrue("(x)->()" in str(f))

    f = ca.Function('f',[],[x],[],["y"], {"allow_free": True})
    self.assertTrue("()->(y)" in str(f))

    f = ca.Function('f',[x],[x**2],["x"],["y"])
    self.assertTrue("(x)->(y)" in str(f))

    f = ca.Function('f',[x,y],[x**2,x*y],["x","y"],["w","z"])
    self.assertTrue("(x,y)->(w,z)" in str(f))

  def test_fold(self):
    x = ca.SX.sym("x",2,2)

    f = ca.Function("f",[x],[ca.sin(x)])

    F = f.fold(10)

    x0 = x
    for i in range(10):
      x = f(x)
    Fref = ca.Function("f",[x0],[x])

    self.checkfunction(F,Fref,inputs=[ca.DM([[1,2],[3,7]])])


  @memory_heavy()
  def test_thread_safety(self):
    x = ca.MX.sym('x')
    y = ca.MX.sym('y')
    f = ca.Function('f', [x, y], [x ** 2 - y])
    for rf in ["newton","fast_newton"]:
      finv = ca.rootfinder('finv', rf, f)

      finv_par = finv.map(50,"unroll").map(4, 'thread',4)
      res = finv_par(numpy.ones(200), numpy.linspace(0, 10, 200))
      self.checkarray(ca.norm_inf(res.T-ca.sqrt(numpy.linspace(0, 10, 200))),0, digits=5)

      finv_par = finv.map(200, 'thread',4)
      res = finv_par(numpy.ones(200), numpy.linspace(0, 10, 200))
      self.checkarray(ca.norm_inf(res.T-ca.sqrt(numpy.linspace(0, 10, 200))),0, digits=5)

  def test_mapped_eval(self):
      x = ca.SX.sym('x')
      y = ca.SX.sym('y', 2)
      f = ca.Function('f', [x,y], [ca.sin(x)*y], ['x','y'], ['r'])
      r1 = f(1, ca.DM([3,4]))
      r2 = f(2, ca.DM([3,4]))
      r3 = f(3, ca.DM([3,4]))
      r_all = f(ca.DM([[1,2,3]]), ca.DM([3,4]))
      self.checkarray(r_all, ca.horzcat(r1,r2,r3), "Mapped evaluation (DM)")

      z = ca.MX.sym("z", 1, 3);
      rz = f(z, ca.DM([3,4]))
      F = ca.Function('F', [z], [rz]);
      r_mx = F(ca.DM([[1,2,3]]))
      self.checkarray(r_all, r_mx, "Mapped evaluation (MX)")

  def test_default_arg(self):
      x = ca.MX.sym("x")
      y = ca.MX.sym("y")
      z = ca.MX.sym("z")

      f = ca.Function('f',[x,y,z],[x,y,z],["x","y","z"],["a","b","c"],{"default_in": [1,2,3]})
      self.check_codegen(f,{"x":5,"z":3},main=True)

  def test_factory_inherit_options(self):
      x = ca.MX.sym("x",5)


      for op in ["grad:f:x","jac:f:x","triu:hess:f:x:x"]:
        f = ca.Function("f",[x],[ca.dot(x,x)],["x"],["f"],{"verbose": True})

        with capture_stdout() as out:
          fgrad = f.factory("fgrad",["x"],[op])

        self.assertTrue("::init" in out[0])


        f = ca.Function("f",[x],[ca.dot(x,x)],["x"],["f"],{"enable_fd":True,"enable_forward":False,"enable_reverse":False,"print_in":True,"fd_method":"central","fd_options":{"h": 1e-7,"h_iter":False}})

        fgrad = f.factory("fgrad",["x"],[op])

        with capture_stdout() as out:
          fgrad(0)

        self.assertTrue("-1e-07," in out[0] or "-1e-007," in out[0] )
        self.assertTrue("1e-07," in out[0] or "1e-007," in out[0] )

  @requires_nlpsol("ipopt")
  @requiresPlugin(ca.Importer,"shell")
  def test_inherit_jit_options(self):

    x = ca.MX.sym("x")

    f = ca.Function('Q',[x],[x**2],{"jit":True,"compiler":"shell","verbose":True})

    [g] = f.call([ca.sin(x)],False,True)
    g = ca.cos(g)
    y = ca.MX.sym("x")
    with capture_stdout() as out:
      solver = ca.nlpsol("solver","ipopt",{"x": x,"f": x**2, "g": g})

    import re
    found = set([m.group(1) for m in re.finditer(r"Compiling function '(\w+?)'", out[0])])
    self.assertTrue(len(found)==1)
    self.assertTrue("fwd1_Q" in found)

    with capture_stdout() as out:
      self.check_serialize(solver,inputs={"x0": 2})

    found = set([m.group(1) for m in re.finditer(r"Compiling function '(\w+?)'", out[0])])
    self.assertTrue(len(found)==2)
    self.assertTrue("Q" in found)
    self.assertTrue("fwd1_Q" in found)
    
  @requiresPlugin(ca.Importer,"shell")
  def test_jit_directory(self):
  
    print(ca.CasadiMeta.feature_list())
  
    print( "ghc-filesystem" not in ca.CasadiMeta.feature_list())

    have_mkstemps = "HAVE_MKSTEMPS" in ca.CasadiMeta.compiler_flags()    
    have_simple_mkstemps = "HAVE_SIMPLE_MKSTEMPS" in ca.CasadiMeta.compiler_flags()
    
    print(have_mkstemps,have_simple_mkstemps)
    create_dirs = "ghc-filesystem" not in ca.CasadiMeta.feature_list()
  
    abs_temp_dir = tempfile.gettempdir()
    
    rel_temp_dir = "jit_dir_test"
    
    
    cwd = os.getcwd()
    
    import shutil
    if os.path.exists(rel_temp_dir):
        shutil.rmtree(rel_temp_dir)
        
    if not os.path.exists("park"):
        os.makedirs("park")
    
    if create_dirs:
        os.makedirs(os.path.join(rel_temp_dir,"foo"))
        if not os.path.exists("foo"):
            os.makedirs("foo")
        if os.name=='nt':
          return # cl.exe not available anyway
        if not os.path.exists("/tmp/foo"):
          os.makedirs("/tmp/foo")
            
    for temp_dir in ["./", rel_temp_dir,abs_temp_dir]:
        ca.GlobalOptions.setTempWorkDir(temp_dir)
        
        
        for directory in ["", "foo", abs_temp_dir]:
            if create_dirs and directory=="foo" and sys.platform in ["win32","darwin"]: continue
                
            print("directory",directory)
            x = ca.MX.sym("x")
            f = None

            if not have_mkstemps and not have_simple_mkstemps:
                with self.assertInException("custom temporary directory"):
                    ca.Function('f',[x],[x**2],{"jit":True,"compiler":"shell", "jit_options": {"verbose":True, "directory": directory}})
                continue
            if True:
                if os.path.isabs(directory):
                    dir = directory
                else:
                    dir = os.path.join(temp_dir,directory)
                    
                gc.collect()
                for e in glob.glob(os.path.join(dir,"tmp_casadi*")):
                    os.remove(e)
                
                f = ca.Function('f',[x],[x**2],{"jit":True,"compiler":"shell", "jit_options": {"verbose":True, "directory": directory}})
                print(dir)
                # All files in the correct location
                self.assertTrue(len(glob.glob(os.path.join(dir,"tmp_casadi*")))>1)
                if "ghc-filesystem" in ca.CasadiMeta.feature_list():
                    os.chdir("park")
                f = None
                gc.collect()
                os.chdir(cwd)
                # All files removed
                self.assertTrue(len(glob.glob(os.path.join(dir,"tmp_casadi*")))==0)
        
        ca.GlobalOptions.setTempWorkDir("./")

  def test_custom_jacobian(self):
    x = ca.MX.sym("x")
    p = ca.MX.sym("p")
    J = ca.Function("jac_Q", [x, p, ca.MX(1,1), ca.MX(1,1)], [x*pi, 1, 1, 1],
        ['x', 'p', 'out_r', 'out_s'], ['jac_r_x', 'jac_r_p', 'jac_s_x', 'jac_s_p'])

    f = ca.Function('Q', [x, p], [x**2, 2*x*p], ['x', 'p'], ['r', 's'],
            dict(custom_jacobian = J, jac_penalty = 0))

    J = None

    [g1,g2] = f.call([x,p],False,True)
    JJ = ca.Function('JJ',[x,p],[ca.jacobian(g1+3*x,x)])

    self.checkarray(JJ(5,2),5*pi+3)

  def test_dump(self):
    x = ca.MX.sym("x",ca.Sparsity.lower(3))
    y = ca.MX.sym("y",0,0)
    z = ca.MX.sym("z",2,2)

    for fmt in ["txt","mtx"]:
      f = ca.Function("f",[x,y,z],[2*x,2*z],["x","y","z"],["a","c"],{"dump":True,"dump_in":True,"dump_out":True,"dump_format":fmt})
      ins = [ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),ca.DM(),ca.DM([[1,3],[4,5]])]
      out = f(*ins)

      F = ca.Function.load("f.casadi")

      x_num = ca.DM.from_file("f.000000.in.x." + fmt)
      y_num = ca.DM.from_file("f.000000.in.y." + fmt)
      z_num = ca.DM.from_file("f.000000.in.z." + fmt)

      self.checkarray(x_num,ins[0])
      self.checkarray(y_num,ins[1])
      self.checkarray(z_num,ins[2])

      a_num = ca.DM.from_file("f.000000.out.a." + fmt)
      c_num = ca.DM.from_file("f.000000.out.c." + fmt)

      self.checkarray(a_num,out[0])
      self.checkarray(c_num,out[1])


      F.generate_in("test_in.txt", [x_num,y_num,z_num])
      F.generate_out("test_out.txt", F.call([x_num,y_num,z_num]))
      X = ca.DM.from_file("f.000000.in.txt")
      A = ca.DM.from_file("f.000000.out.txt")

      Xr = ca.DM.from_file("test_in.txt")
      Ar = ca.DM.from_file("test_out.txt")

      self.checkarray(Xr,X)
      self.checkarray(Ar,A)

    if args.run_slow and "ghc-filesystem" in ca.CasadiMeta.feature_list():
      import shutil
      for d in ["dump_fn", "dump_fn_codegen"]:
        if os.path.exists(d):
          shutil.rmtree(d)

      # Codegen test for dump_in/dump_out
      x = ca.MX.sym("x", ca.Sparsity.lower(3))
      z = ca.MX.sym("z", 2, 2)
      f = ca.Function("f", [x, z], [2*x, 2*z], ["x", "z"], ["a", "c"],
                   {"dump_in": True, "dump_out": True, "dump_dir": "dump_fn"})

      ins = [ca.sparsify(ca.DM([[1, 0, 0], [2, 4, 0], [7, 8, 9]])), ca.DM([[1, 3], [4, 5]])]
      f(*ins)

      self.check_codegen(f, inputs=ins, std="c99", main=True,
                         opts={"dump_dir_suffix": "_codegen"})

      for name in ["f.000000.in.txt", "f.000000.in.x.mtx", "f.000000.in.z.mtx",
                    "f.000000.out.a.mtx", "f.000000.out.c.mtx", "f.000000.out.txt"]:
        self.file_equal(os.path.join("dump_fn", name),
                        os.path.join("dump_fn_codegen", name))

      for d in ["dump_fn", "dump_fn_codegen"]:
        if os.path.exists(d):
          shutil.rmtree(d)

  def test_print(self):
    import re
    x = ca.MX.sym("x",ca.Sparsity.lower(3))
    z = ca.MX.sym("z",2,2)
    f = ca.Function("f",[x,z],[2*x,2*z],["x","z"],["a","c"],
                 {"print_in":True,"print_out":True,"print_canonical":True})

    ins = [ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),ca.DM([[1,3],[4,5]])]

    with capture_stdout() as out:
      result = f(*ins)
    vm_output = out[0]

    self.checkarray(result[0],2*ins[0])
    self.checkarray(result[1],2*ins[1])

    # Strip pointer addresses from VM output for comparison
    vm_lines = [re.sub(r' \(0x[0-9a-f]+\)','',line)
                for line in vm_output.strip().split('\n')]
    self.assertTrue(len(vm_lines)>=6)

    self.check_codegen(f,inputs=ins,std="c99",main=True,main_output_check=False)
    if args.run_slow:
      with open(f.name()+"_out.txt","r") as fh:
        cg_output = fh.read()
      # Extract print lines (skip main output values line)
      cg_lines = [line for line in cg_output.strip().split('\n')
                  if line.startswith("Function") or line.startswith("Input")
                     or line.startswith("Output")]
      self.assertEqual(vm_lines,cg_lines)

  def test_eval_shapes(self):
    x = ca.MX.sym("x",ca.Sparsity.lower(3))
    y = ca.MX.sym("y",3,1)
    z = ca.MX.sym("z",2,2)
    f = ca.Function("f",[x,y,z],[2*x,2*y,2*z],{"default_in":[1,2,3]})

    for ins, ref, ref_flat in [
      #normal
      ([ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),ca.DM([9,8,7]),ca.DM([[1,3],[4,5]])],
       [2*ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])), 2*ca.DM([9,8,7]), 2*ca.DM([[1,3],[4,5]])],
       [1,2,7,4,8,9,9,8,7,1,4,3,5]),

      #Project
      ([ca.sparsify(ca.DM([[1,0,0],[2,0,0],[7,8,9]])),ca.DM([9,8,7]),ca.DM([[1,3],[4,5]])],
       [2*ca.sparsify(ca.DM([[1,0,0],[2,1e-30,0],[7,8,9]])),2*ca.DM([9,8,7]),2*ca.DM([[1,3],[4,5]])],
       [1,2,7,0,8,9,9,8,7,1,4,3,5]),

      #Project away
      ([ca.sparsify(ca.DM([[1,0,1],[2,4,0],[7,8,9]])),ca.DM([9,8,7]),ca.DM([[1,3],[4,5]])],
       [2*ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])), 2*ca.DM([9,8,7]), 2*ca.DM([[1,3],[4,5]])],
       [1,2,7,4,8,9,9,8,7,1,4,3,5]),

      #Scalar expansion
      ([ca.DM(4),5,7],
       [2*ca.sparsify(ca.DM([[4,0,0],[4,4,0],[4,4,4]])), 2*ca.DM([5,5,5]), 2*ca.DM([[7,7],[7,7]])],
       [4,4,4,4,4,4,5,5,5,7,7,7,7]),


      #Repmat
      ([ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),ca.DM([9,8,7]),ca.DM([1,3])],
       [2*ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])), 2*ca.DM([9,8,7]), 2*ca.DM([[1,1],[3,3]])],
       [1,2,7,4,8,9,9,8,7,1,3,1,3]),

      #npar
      ([ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),ca.DM([[9,2],[8,4],[7,4]]),ca.DM([[1,3],[4,5]])],
       [2*ca.sparsify(ca.DM([[1,0,0,1,0,0],[2,4,0,2,4,0],[7,8,9,7,8,9]])), 2*ca.DM([[9,2],[8,4],[7,4]]), 2*ca.DM([[1,3,1,3],[4,5,4,5]])],
       None),


      # growing npar
      ([ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),ca.DM([[9,2],[8,4],[7,4]]),ca.DM([[1,3,4,5,1,1,1,1],[4,5,7,8,3,3,3,3]])],
       [2*ca.repmat(ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),1,4), 2*ca.repmat(ca.DM([[9,2],[8,4],[7,4]]),1,2), 2*ca.DM([[1,3,4,5,1,1,1,1],[4,5,7,8,3,3,3,3]])],
       None),

      # growing npar
      ([ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),ca.DM([[9,2,3,4],[8,4,7,8],[7,4,9,10]]),ca.DM([[1,3,4,5],[4,5,7,8]])],
       None,
       None),

      #Transpose
      ([ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),ca.DM([9,8,7]).T,ca.DM([[1,3],[4,5]])],
       [2*ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),2*ca.DM([9,8,7]), 2*ca.DM([[1,3],[4,5]])],
       [1,2,7,4,8,9,9,8,7,1,4,3,5]),

      #Null: note default_in applies only to dict style
      ([ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),ca.DM(),ca.DM([1,3])],
       [2*ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),2*ca.DM([0,0,0]), 2*ca.DM([[1,1],[3,3]])],
       [1,2,7,4,8,9,0,0,0,1,3,1,3]),

      #Null: note default_in applies only to dict style
      ([ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),ca.DM(3,0),ca.DM([1,3])],
       [2*ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),2*ca.DM([0,0,0]), 2*ca.DM([[1,1],[3,3]])],
       [1,2,7,4,8,9,0,0,0,1,3,1,3]),

      #Null: note default_in applies only to dict style
      ([ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),ca.DM(0,3),ca.DM([1,3])],
       [2*ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),2*ca.DM([0,0,0]), 2*ca.DM([[1,1],[3,3]])],
       [1,2,7,4,8,9,0,0,0,1,3,1,3]),

      #Null: should be forbidden
      ([ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),ca.DM(0,5),ca.DM([1,3])],
       [2*ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),2*ca.DM([0,0,0]), 2*ca.DM([[1,1],[3,3]])],
       [1,2,7,4,8,9,0,0,0,1,3,1,3]),

      #Null: should be forbidden
      ([ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),ca.DM(8,0),ca.DM([1,3])],
       [2*ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),2*ca.DM([0,0,0]), 2*ca.DM([[1,1],[3,3]])],
       [1,2,7,4,8,9,0,0,0,1,3,1,3]),



      # Wrong
      ([ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),ca.DM([9,8,7]).T,ca.DM([1,3]).T],None,None),
      ([ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9],[4,5,6]])),ca.DM([9,8,7]),ca.DM([1,3])],None,None),
      ([ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),ca.DM([[9,2,1],[8,4,1],[7,4,1]]),ca.DM([[1,1,1,3],[3,4,5,6]])],None,None),
      ([ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),ca.DM([9,8,7]),ca.DM([[1,1,1],[3,4,5]])],None,None),
      ([ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])),ca.DM([9,8,7,6]),ca.DM([[1,3],[4,5]])],None,None),

      ]:

      if ref is None:
        with self.assertInException("mismatching shape"):
          res = f.call(ins)
      else:
        res = f.call(ins)

        for i in range(f.n_in()):
          self.assertTrue(res[i].sparsity()==ref[i].sparsity())
          self.checkarray(res[i],ref[i])

        if ref_flat is None:
          with self.assertInException("mismatching shape"):
            f.nz_from_out(res)
        else:
          res_flat = f.nz_from_out(res)

          self.checkarray(res_flat,[i*2 for i in ref_flat])

          res_re = f.nz_to_out(res_flat)

          for i in range(f.n_out()):
            self.assertTrue(res_re[i].sparsity()==ref[i].sparsity())
            self.checkarray(res_re[i],ref[i])

          f.generate_out("test.txt",res)
          res2 = f.generate_in("test.txt")

          for i in range(f.n_out()):
            self.assertTrue(res2[i].sparsity()==ref[i].sparsity())
            self.checkarray(res2[i],ref[i])

      if ref_flat is None:
        with self.assertInException("mismatching shape"):
          res = f.nz_from_in(ins)

        with self.assertInException("mismatching shape"):
          f.generate_in("test.txt",ins)
      else:
        assert ref is not None  # by construction: ref_flat is not None implies ref is not None
        res = f.nz_from_in(ins)
        self.checkarray(res,ref_flat)

        in_re = f.nz_to_in(res)
        for i in range(f.n_in()):
          self.assertTrue(in_re[i].sparsity()==ref[i].sparsity())
          self.checkarray(in_re[i],ref[i]/2)

        f.generate_in("test.txt",ins)
        ins2 = f.generate_in("test.txt")

        for i in range(f.n_in()):
          self.assertTrue(ins2[i].sparsity()==ref[i].sparsity())
          self.checkarray(ins2[i],ref[i]/2)


    #Null dicts
    ins = {"i0": ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])), "i2": ca.DM([1,3])}
    res = f.call(ins)

    self.checkarray(res["o0"],2*ca.sparsify(ca.DM([[1,0,0],[2,4,0],[7,8,9]])))
    self.checkarray(res["o1"],2*ca.DM([2,2,2]))
    self.checkarray(res["o2"],2*ca.DM([[1,1],[3,3]]))

    res = f.nz_from_in(f.convert_in(ins))
    self.checkarray(res,[1,2,7,4,8,9,2,2,2,1,3,1,3])

  def test_issue3977(self):
    # A 0-by-N argument to a 0-by-M input must be recognised as a parallel/repmat
    # call (N a multiple of M), not collapsed to the input shape as a generic empty.
    for X in [ca.SX, ca.MX]:
      x = X.sym("x", 0, 1)
      f = ca.Function("f", [x], [x**2])
      self.assertEqual(f(ca.DM(0, 3)).shape, (0, 3))
      self.assertEqual(f(ca.DM(0, 1)).shape, (0, 1))

      x = X.sym("x", 0, 2)
      f = ca.Function("f", [x], [x**2])
      self.assertEqual(f(ca.DM(0, 6)).shape, (0, 6))

    # Empty argument against a non-degenerate input still acts as a zero placeholder.
    x = ca.SX.sym("x", 3, 1)
    f = ca.Function("f", [x], [x**2])
    self.assertEqual(f(ca.DM(0, 0)).shape, (3, 1))
    self.assertEqual(f(ca.DM(0, 1)).shape, (3, 1))

  def test_convert_in(self):
    x = ca.MX.sym("x")
    y = ca.MX.sym("y")
    z = ca.MX.sym("z")
    f = ca.Function("f",[x,y,z],[2*x,2*y,2*z,x+y],["x","y","z"],["a","b","c","d"],{"default_in": [1,2,3]})

    res = f.convert_in({"x":4})
    self.checkarray(res[0],4)
    self.checkarray(res[1],2)
    self.checkarray(res[2],3)

    with self.assertInException("xy"):
      f.convert_in({"xy":4})

    res = f.convert_in([7,3,4])
    self.checkarray(res["x"],7)
    self.checkarray(res["y"],3)
    self.checkarray(res["z"],4)

    with self.assertInException("Incorrect"):
      f.convert_in([7,3])

    res = f.convert_out({"b":4})
    self.checkarray(res[0],np.nan)
    self.checkarray(res[1],4)
    self.checkarray(res[2],np.nan)
    self.checkarray(res[3],np.nan)

    with self.assertInException("xy"):
      f.convert_out({"xy":4})

    res = f.convert_out([1,2,3,4])
    self.checkarray(res["a"],1)
    self.checkarray(res["b"],2)
    self.checkarray(res["c"],3)
    self.checkarray(res["d"],4)

    with self.assertInException("Incorrect"):
      f.convert_out([1,2,3])

  @requires_nlpsol("ipopt")
  def test_ad_weight_sp(self):
    x = ca.MX.sym("x",4)

    f = ca.Function('f',[x],[x])

    with self.assertOutput("1 forward sweeps needed","1 reverse sweeps needed"):
      ca.jacobian(f(ca.sin(x)),x,{"helper_options": {"verbose":True, "ad_weight_sp":0}})
    with self.assertOutput("1 reverse sweeps needed","1 forward sweeps needed"):
      ca.jacobian(f(ca.sin(x)),x,{"helper_options": {"verbose":True, "ad_weight_sp":1}})


    x = ca.MX.sym("x",4)

    nlp = {"x":x,"f":ca.dot(x,x),"g": ca.sin(x)}

    with self.assertOutput("1 forward sweeps needed","1 reverse sweeps needed"):
      solver = ca.nlpsol("solver","ipopt",nlp,{"common_options": {"helper_options": {"verbose":True, "ad_weight_sp":0}}})
    with self.assertOutput("1 reverse sweeps needed","1 forward sweeps needed"):
      solver = ca.nlpsol("solver","ipopt",nlp,{"common_options": {"helper_options": {"verbose":True, "ad_weight_sp":1}}})


    f = ca.Function('f',[x],[x],{"verbose":True,"ad_weight_sp":1})
    with self.assertOutput("1 reverse sweeps needed","1 forward sweeps needed"):
      ca.jacobian(f(ca.sin(x)),x)

    x = ca.MX.sym("x",100)

    f = ca.Function('f',[x],[x[0]],{"verbose":True,"ad_weight_sp":0})

    with self.assertOutput("2 forward sweeps needed","1 reverse sweeps needed"):
      ca.jacobian(f(ca.sin(x)),x)

  def test_map_exception(self):
    x = ca.MX.sym("x",4)
    y = ca.MX.sym("y",4)
    f = ca.Function("f",[x],[x+y],{"allow_free":True})

    if "CASADI_WITH_THREAD" in ca.CasadiMeta.compiler_flags():
      message = "Evaluation failed"
    else:
      message = "since variables [y] are free"

    with self.assertInException(message):
      F = f.map(4,"thread",2)
      F(3)
      
  def test_DM_arg(self):
    f = ca.Function('f',[ca.DM(0,1)],[])
    print(f)
    f = ca.Function('f',[ca.DM(0,1),ca.SX.sym("x")],[])
    print(f)
    f = ca.Function('f',[ca.DM(0,1),ca.MX.sym("x")],[])
    print(f)
    self.assertTrue(isinstance(ca.vertcat(),ca.DM))
    self.assertEqual(ca.vertcat().shape,(0,1))
    self.assertTrue(isinstance(ca.horzcat(),ca.DM))
    self.assertEqual(ca.horzcat().shape,(1,0))
    self.assertTrue(isinstance(ca.veccat(),ca.DM))
    self.assertEqual(ca.veccat().shape,(0,1))
    self.assertTrue(isinstance(ca.diagcat(),ca.DM))
    self.assertEqual(ca.diagcat().shape,(0,0))
    self.assertTrue(isinstance(ca.vcat([]),ca.DM))
    self.assertEqual(ca.vcat([]).shape,(0,1))
    self.assertTrue(isinstance(ca.hcat([]),ca.DM))
    self.assertEqual(ca.hcat([]).shape,(1,0))
    self.assertTrue(isinstance(ca.vvcat([]),ca.DM))
    self.assertEqual(ca.vvcat([]).shape,(0,1))
    self.assertTrue(isinstance(ca.dcat([]),ca.DM))
    self.assertEqual(ca.dcat([]).shape,(0,0))
    
    
  def test_nondiff(self):

    for X in [ca.SX,ca.MX]:

      x = X.sym("x",2,2)
      y = X.sym("y",2,2)

      xn = ca.DM([[1,3],[0.1,2]])
      yn = ca.DM([[4,5],[6,7]])

      options = {"is_diff_in":[True,False],"is_diff_out":[True,False]}
      f = ca.Function("f",[x,y],[ca.sin(x+y),x*y],options)
      
      if X is ca.MX:
        fexpanded = f.expand()
        r1 = str(f).split(" ")[0]
        r2 = str(fexpanded).split(" ")[0]
        self.assertEqual(r1,r2)



      F = ca.Function("F",[x,y],f(ca.cos(x),(x*y)**2),["x","y"],["z","zz"],options)


      for ff in [F.forward(1).forward(1),F.forward(1).reverse(1),F.reverse(1).forward(1),F.reverse(1).reverse(1)]:
        s = str(ff)
        self.assertTrue("y" not in s.replace("_y[2x2,0nz]","foo")[len("fwd1_adj1_F:(x[2x2],y[2x2]"):])


      ye = X(2,2);

      G = ca.Function("G",[x,ye],[f(ca.cos(x),(x*yn)**2)[0], f(ca.cos(xn),(xn*yn)**2)[1]],["x","y"],["z","zz"])
      self.checkfunction(F,G,inputs=[xn,yn],evals=1)

      for f1 in [lambda f: f.forward(1), lambda f: f.reverse(1)]:
         Gf = f1(G)
         Ff = f1(F)

         arg = [xn,yn]+[ca.DM.rand(Gf.sparsity_in(i)) for i in range(Gf.n_in())][2:]
         for a,b in zip(Gf.call(arg),Ff.call(arg)):
            self.checkarray(a,b)

      for f1 in [lambda f: f.forward(1), lambda f: f.reverse(1)]:
        for f2 in [lambda f: f.forward(1), lambda f: f.reverse(1)]:
           Gf = f1(f2(G))
           Ff = f1(f2(F))

           arg = [xn,yn]+[ca.DM.rand(Gf.sparsity_in(i)) for i in range(Gf.n_in())][2:]
           for a,b in zip(Gf.call(arg),Ff.call(arg)):
             self.checkarray(a,b)

  def test_non_diff_call(self):
    x = ca.MX.sym("x",5)
    y = ca.MX.sym("x",5,5)

    f = ca.Function("f",[x,y],[y @ x,y.T @ x],{"never_inline":True,"is_diff_in":[False,True],"is_diff_out":[False,True]})
    
    

    x = ca.SX.sym("x",5)
    y = ca.SX.sym("x",5,5)

    g = ca.Function('g',[x,y],list(f(x,y)))


    X = ca.MX.sym("x",5)
    Y = ca.MX.sym("x",5,5)

    G = ca.Function('g',[x,y],list(f(x,y)))
    
    ca.DM.rng(1)
    x0=ca.DM.rand(5)
    y0=ca.DM.rand(5,5)
    
    self.checkfunction(g,G,inputs=[x0,y0])
    

  @memory_heavy()
  def test_inline_linear_interpolant(self):

    do_inline = False

    np.random.seed(0)

    N = 3
    M = 4

    d_knots = [list(np.linspace(0,1,N)),list(np.linspace(0,1,M))]

    data0 = np.random.random([len(e) for e in d_knots])
    data1 = np.random.random([len(e) for e in d_knots])
    r = np.meshgrid(*d_knots,indexing='ij')

    xyz = np.vstack(list(e.ravel(order='F') for e in r)).ravel(order='F')

    d_flat0 = data0.ravel(order='F')
    d_flat1 = data1.ravel(order='F')

    data = np.vstack((data0.ravel(order='F'),data1.ravel(order='F'))).ravel(order='F')
    d_flat = data.ravel(order='F')


    LUT_param = ca.interpolant('name','linear',[N,M],2, {"lookup_mode": ["exact","linear"]})

    only_fd = {"enable_fd": True, "enable_forward": False, "enable_reverse": False}
    options = {"ad_weight_sp":-1}  # type: dict
    options.update(only_fd)
    options["forward_options"] = only_fd
    options["reverse_options"] = only_fd
    #options["jacobian_options"] = only_fd

    print(options)
    LUT_param_ref = LUT_param.wrap_as_needed(options)
    J_ref = LUT_param_ref.jacobian()

    d_knots_cat = ca.vcat(d_knots[0]+d_knots[1])

    inputs = [[0.2,0.333],d_knots_cat,d_flat]

    LUT_param = ca.interpolant('name','linear',[N,M],2,{"inline": True, "lookup_mode": ["exact","linear"]})
    J = LUT_param.jacobian()

    self.checkarray(LUT_param(*inputs),LUT_param_ref(*inputs))
    inputs+= [0]
    print("J_ref",J_ref,J_ref(*inputs))
    self.checkfunction(J,J_ref,inputs=inputs,digits=8,digits_sens=1,evals=1)


  def test_issue3079(self):  
    ca.interp1d([1,2,3], ca.MX([4,5,6]), [1.1], 'linear')

  def test_functionbuffer(self):
    A_ = np.random.random((4,4))
    B_ = np.zeros((4,4))

    A = ca.MX.sym("A",4,4)

    f = ca.Function("f",[A],[A*10])


    [buf,f_eval] = f.buffer()
    buf.set_arg(0, memoryview(A_))  # pyright: ignore[reportCallIssue,reportArgumentType]
    buf.set_res(0, memoryview(B_))  # pyright: ignore[reportCallIssue,reportArgumentType]

    f_eval()

    self.checkarray(B_, 10*A_)

    self.assertEqual(buf.ret(), 0)

  @requires_conic("osqp")
  @requiresPlugin(ca.Importer,"shell")
  def test_jit_buffer_eval(self):
    opti = ca.Opti('conic')

    x = opti.variable()

    opti.minimize((x-3)**2)

    opti.subject_to(x>=0)

    opti.subject_to(x<=5)


    opti.solver('osqp')


    f = opti.to_function('f',[],[x])


    f.generate('f.c')

    libdir = ca.GlobalOptions.getCasadiPath()
    includedir = ca.GlobalOptions.getCasadiIncludePath()

    if os.name=='nt':
      f = opti.to_function('f',[],[x],{"jit":True,"compiler":"shell","jit_options":{"compiler_flags":["/I"+includedir],"linker_flags":["/LIBPATH:"+libdir,"osqp.lib"],"verbose":True}, "print_time": True})
    else:
      f = opti.to_function('f',[],[x],{"jit":True,"compiler":"shell","jit_options":{"compiler_flags":["-I"+includedir,"-L"+libdir,"-losqp"],"linker_flags":["-I"+includedir,"-L"+libdir,"-losqp"],"verbose":True}, "print_time": True})
    [buf,trigger] = f.buffer()

    a = np.array([1.0])
    buf.set_res(0, memoryview(a))  # pyright: ignore[reportCallIssue,reportArgumentType]

    trigger()

    # buffer jit eval bypasses timings
    with self.assertInException("No stats available"):
        f.stats()
    self.checkarray(a,3)

  def test_codegen_inf_nan(self):
    x = ca.MX.sym("x")
    f = ca.Function("F",[x],[x+inf])
    self.check_codegen(f,inputs=[1],std="c99",main=True)

    x = ca.MX.sym("x")
    f = ca.Function("F",[x],[x+np.nan])
    self.check_codegen(f,inputs=[1],std="c99",main=True)

    x = ca.MX.sym("x")
    f = ca.Function("F",[x],[x+ca.vertcat(inf,np.nan,-inf)])
    self.check_codegen(f,inputs=[1],std="c99",main=True)

  def test_codegen_with_mem(self):
    x = ca.MX.sym("x")
    f = ca.Function("F",[x],[3*x])
    self.check_codegen(f,inputs=[1],main=True,opts={"with_mem":True,"force_canonical":True},definitions=["inline=\"\""])
    self.check_codegen(f,inputs=[1],main=True,opts={"with_mem":True,"with_header":True,"force_canonical":True},definitions=["inline=\"\""])

  def test_codegen_scalars_bug(self):
    x = ca.MX.sym("x")
    z = 3*x/ca.sin(x)
    f = ca.Function("F",[x],[z,z/x],{"live_variables":False})
    self.check_codegen(f,inputs=[1],opts={"codegen_scalars":True})

  def test_bug_codegen_logical(self):
    a = ca.MX([1,0,0])
    b = ca.MX([1,1,0])
    c = ca.logic_or(a,b)
    f = ca.Function("f",[],[c])
    self.check_codegen(f,inputs=[])

  def test_jit_serialize(self):
    if not args.run_slow: return
    if sys.platform=="darwin": return

    def test_cases(jit_serialize):
      opts = {"jit":True, "compiler": "shell", "jit_options": {"verbose":True}, "verbose":False, "jit_serialize": jit_serialize}
      x = ca.MX.sym("x")
      yield lambda : ca.Function('f',[x],[(x-3)**2],opts)


      x = ca.MX.sym("x", 2)
      f = x[0]*x[0] + x[1]*x[1]

      if ca.has_nlpsol("ipopt"):
        yield lambda : ca.nlpsol("solver", "ipopt", {"x": x, "f": f},opts)



    for case in test_cases("source"):

      with self.assertOutput(["jit_tmp"],[]):
        f = case()

      f.save('f.casadi')

      with self.assertOutput(["jit_tmp"],[]):
        g = ca.Function.load('f.casadi')

      self.checkfunction_light(f, g, inputs={})

      g.save('f.casadi')
      with self.assertOutput(["jit_tmp"],[]):
        g = ca.Function.load('f.casadi')

      self.checkfunction_light(f, g, inputs={})



    for case in test_cases("link"):

      with self.assertOutput(["jit_tmp"],[]):
        f = case()


      f.save('f.casadi')

      with self.assertOutput([],["jit_tmp"]):
        g = ca.Function.load('f.casadi')

      self.checkfunction_light(f, g, inputs={})

      g.save('f.casadi')
      with self.assertOutput([],["jit_tmp"]):
        g = ca.Function.load('f.casadi')

      self.checkfunction_light(f, g, inputs={})

      f = None
      g = None
      if os.name!='nt': # Workaround for known bug #3039
        with self.assertInException("No such file"):
          g = ca.Function.load('f.casadi')


    for case in test_cases("embed"):

      with self.assertOutput(["jit_tmp"],[]):
        f = case()

      f.save('f.casadi')
      with self.assertOutput([],["jit_tmp"]):
        g = ca.Function.load('f.casadi')
      g()
      f = None
      with self.assertOutput([],["jit_tmp"]):
        g = ca.Function.load('f.casadi')

      g()
      g.save('f.casadi')
      with self.assertOutput([],["jit_tmp"]):
        g = ca.Function.load('f.casadi')
      g()
      g = None
      with self.assertOutput([],["jit_tmp"]):
        g = ca.Function.load('f.casadi')


  def test_map_get_function(self):
    x = ca.MX.sym("x")

    g = ca.Function("g",[x],[x**2])
    ff = g.map(5)

    self.assertTrue(len(ff.get_function())==1)
    f2 = ff.get_function("f")

    self.checkfunction_light(g, f2, inputs=[3])

  def test_find_function(self):
    for X in [ca.MX,ca.SX]:
        x = X.sym("x")
        
        """
           k
           j
         i    g
        h  g      
        """

        g = ca.Function("g",[x],[x**2],{"never_inline":True})
        
        h = ca.Function("h",[x],[g(x**2)],{"never_inline":True})
        
        i = ca.Function("i",[x],[h(x**2)+g(x**2)],{"never_inline":True})
        
        j = ca.Function("j",[x],[i(x**2)+g(x**2)],{"never_inline":True})
        
        k = ca.Function("k",[x],[j(x**2)],{"never_inline":True})
        
        k.generate('f.c')
         
        ref = [j]
        res = k.find_functions(0)
        
        self.assertEqual(len(ref),len(res))
        for a,b in zip(ref,res):  assert(a.name() == b.name())
          
        ref = [j,i,g]
        res = k.find_functions(1)
        
        self.assertEqual(len(ref),len(res))
        for a,b in zip(ref,res):  assert(a.name() == b.name())
        
        ref = [j,i,h,g]
        res = k.find_functions(2)
        
        self.assertEqual(len(ref),len(res))
        for a,b in zip(ref,res):  assert(a.name() == b.name())    

        ref = [j,i,h,g]
        res = k.find_functions()
        
        self.assertEqual(len(ref),len(res))
        for a,b in zip(ref,res):  assert(a.name() == b.name())    


        intg = ca.integrator("intg","rk",{"x":x,"ode":k(x)})
        ref = [intg.get_function("dae"),k,j,i,h,g,intg.get_function("daeF"),intg.get_function("step")]
        res = intg.find_functions()
        print(res)
        
        self.assertEqual(len(ref),len(res))
        for a,b in zip(ref,res):  assert(a.name()==b.name())    

        K = k.map(5)
        
        ref = [k,j,i,h,g]
        res = K.find_functions()
        
        self.assertEqual(len(ref),len(res))
        for a,b in zip(ref,res):  assert(a.name() == b.name())    

  def test_post_expand(self):

    x = ca.MX.sym("x")

    f = ca.Function("f",[x],[x**2])

    print(f)

    f = ca.Function("f",[x],[x**2],{"post_expand":True})
    print(f)
    self.assertTrue(f.is_a("SXFunction"))

    f = ca.Function("f",[x],[x**2],{"post_expand":True,"post_expand_options":{"print_in":True}})
    print(f)
    self.assertTrue(f.is_a("SXFunction"))

    # Check that the option came through
    with self.assertOutput(["Input 0 (i0): 3"],[]):
      f(3)
  
  @requires_conic('osqp')
  def test_memful_main(self):
    c = ca.conic("conic","osqp",{"a":ca.Sparsity.dense(1,2),"h":ca.Sparsity.dense(2,2)},{"print_problem":True,"osqp.verbose":False})
    inputs = {"h": ca.DM([[1,0.2],[0.2,1]]),"g":ca.vertcat(1,2),"a":ca.horzcat(1,1),"lba":-1,"uba":1}
    extra_options = ["-Wno-endif-labels","-Wno-unused-variable","-Wno-strict-prototypes"]
    # No extra options on windows
    if os.name == 'nt':
      extra_options = []
    
    self.check_codegen(c,inputs=c.convert_in(inputs),main=True,valgrind=True,std="c99",extra_options=extra_options,extralibs=["osqp"])
    
    c.generate('f.c')
    
    n_matches = 0
    
    with open("f.c","r") as inp:
        for line in inp.readlines():
            if re.search(r"static int casadi_(\w+)_ref_counter",line):
                n_matches = n_matches +1
    
    self.assertEqual(n_matches,1)
    
    # Check wrapper Function
    for X in [ca.MX,ca.SX]:
    
        g = X.sym("g",2)
        inputs["g"] = g
        
        f = ca.Function('wrapper',[g],[c(**inputs)['x']])
        print(f)
        self.check_codegen(f,inputs=[ca.vertcat(1,2)],main=True,valgrind=True,std="c99",extra_options=extra_options,extralibs=["osqp"])
        f.generate('f.c')
        
        n_matches = 0
        
        with open("f.c","r") as inp:
            for line in inp.readlines():
                if re.search(r"static int casadi_(\w+)_ref_counter",line):
                    n_matches = n_matches +1
        
        # wrapper function should not have it's own reference counting
        self.assertEqual(n_matches,1)

  @memory_heavy()
  def test_codegen_thread_safe(self):
  
    for variant in [0,1]:
        if variant==0 and ca.has_conic("osqp"):
            c = ca.conic("conic","osqp",{"a":ca.Sparsity.dense(1,2),"h":ca.Sparsity.dense(2,2)},{"osqp.verbose":False,"osqp.alpha":1,"osqp.eps_abs":1e-8,"osqp.eps_rel":1e-8,"print_time":False})
            inputs = {"h": ca.DM([[1,0.2],[0.2,1]]),"g":ca.vertcat(1,2),"a":ca.horzcat(1,1),"lba":-1,"uba":1}
            extralibs = ["osqp"]
            symbol = "g"
        elif ca.has_nlpsol("fatrop") and os.name == 'posix':
            x = ca.MX.sym('x',2)
            solver = ca.nlpsol("solver","fatrop",{"x":x,"f":ca.sumsqr(x-10),"g":x[0]+x[1]},{"fatrop.print_level":0,"print_time":False})
            inputs = {"x0": ca.vertcat(1,2),"lbx":ca.vertcat(-5,-4),"ubx":ca.vertcat(5,4),"lbg":-2,"ubg":2}
            extralibs = ["fatrop","blasfeo","m"]
            c= solver
            symbol = "ubx"
        else:
            continue
        extra_options = ["-Wno-endif-labels","-Wno-unused-variable","-Wno-strict-prototypes"]
        # No extra options on windows
        if os.name == 'nt':
          extra_options = []
          
        def get_functions():
            yield (c,inputs)
            # Check wrapper Function
            for X in [ca.MX,ca.SX]:
                local_inputs = dict(inputs)
            
                g = X.sym(symbol,2)
                local_inputs[symbol] = g
                
                yield (ca.Function('wrapper',[g],[c(**local_inputs)['x']],[symbol],["res"]),{symbol: inputs[symbol]})
                
        for f,local_inputs in get_functions():
            for map_type in ["serial","thread","openmp"]:
                F = f.map(2,map_type)
                


                wide_inputs = {}
                for k,v in local_inputs.items():
                    wide_inputs[k] = ca.horzcat(local_inputs[k],local_inputs[k]*1.2)

                if map_type == "openmp":
                    # Set up OpenMP compilation flags
                    if os.name == 'nt':
                        openmp_flags = ["/openmp"]
                    elif sys.platform == 'darwin':
                        # macOS: use -Xpreprocessor -fopenmp and environment variables
                        openmp_flags = ["-Xpreprocessor", "-fopenmp"]
                        # Add CPPFLAGS and LDFLAGS from environment if available (set by CI)
                        cppflags = os.environ.get('CPPFLAGS', '')
                        ldflags = os.environ.get('LDFLAGS', '')
                        if cppflags:
                            openmp_flags.extend(cppflags.split())
                        if ldflags:
                            openmp_flags.extend(ldflags.split())
                        openmp_flags.append("-lomp")
                        openmp_flags.append("-Wno-pedantic") # workaround /opt/homebrew/opt/libomp/include/omp.h:60:9: error: ISO C restricts enumerator values to range of 'int' (2147483648 is too large) [-Werror,-Wpedantic]
                    else:
                        openmp_flags = ["-fopenmp"]

                    self.check_codegen(F,inputs=F.convert_in(wide_inputs),main=True,valgrind=True,std="c99",extra_options=extra_options+openmp_flags,extralibs=extralibs,opts={"thread_safe":True},definitions=["CASADI_MAX_NUM_THREADS=2","CASADI_THREAD_TYPE=CASADI_THREAD_TYPE_OMP"],with_external=False,digits=12)
                elif map_type == "thread":
                    if os.name == "posix":
                        for CASADI_MUTEX_USE_STATIC_INIT in ["0","1"]:
                            self.check_codegen(F,inputs=F.convert_in(wide_inputs),main=True,valgrind=True,std="c99",extra_options=extra_options+["-pthread"],extralibs=extralibs,opts={"thread_safe":True},definitions=["CASADI_MAX_NUM_THREADS=2","CASADI_THREAD_TYPE=CASADI_THREAD_TYPE_POSIX","CASADI_MUTEX_USE_STATIC_INIT=" + CASADI_MUTEX_USE_STATIC_INIT],with_external=False,helgrind=True,digits=12)
                        if sys.platform != "darwin":
                            self.check_codegen(F,inputs=F.convert_in(wide_inputs),main=True,valgrind=True,std="c11",extra_options=extra_options,extralibs=extralibs,opts={"thread_safe":True},definitions=["CASADI_MAX_NUM_THREADS=2","CASADI_THREAD_TYPE=CASADI_THREAD_TYPE_C11"],with_external=False,helgrind=True,digits=12)
                    elif os.name == 'nt':
                        #self.check_codegen(F,inputs=F.convert_in(wide_inputs),main=True,std="c11",extra_options=extra_options,extralibs=extralibs,opts={"thread_safe":True},definitions=["CASADI_MAX_NUM_THREADS=2","CASADI_THREAD_TYPE=CASADI_THREAD_TYPE_C11"],with_external=False,debug_mode=True)
                        for CASADI_MUTEX_USE_STATIC_INIT in ["0","1"]: # Seems to only work when casadi is compiled with visual studio
                            self.check_codegen(F,inputs=F.convert_in(wide_inputs),main=True,std="c99",extra_options=extra_options,extralibs=extralibs,opts={"thread_safe":True},definitions=["CASADI_MAX_NUM_THREADS=2","CASADI_THREAD_TYPE=CASADI_THREAD_TYPE_WINDOWS","CASADI_MUTEX_USE_STATIC_INIT="+CASADI_MUTEX_USE_STATIC_INIT],with_external=False,debug_mode=True,digits=12)
                else:
                    self.check_codegen(f,inputs=f.convert_in(local_inputs),main=True,valgrind=True,std="c99",extra_options=extra_options,extralibs=extralibs,with_external=False,digits=12)
                    self.check_codegen(f,inputs=f.convert_in(dict((k,v*1.2) for k,v in local_inputs.items())),main=True,valgrind=True,std="c99",extra_options=extra_options,extralibs=extralibs,with_external=False,digits=12)
                    self.check_codegen(F,inputs=f.convert_in(wide_inputs),main=True,valgrind=True,std="c99",extra_options=extra_options,extralibs=extralibs,with_external=False,digits=12)
            

  @requires_conic('osqp')
  def test_memful_external(self):
    if not args.run_slow: return
    c = ca.conic("conic","osqp",{"a":ca.Sparsity.dense(1,2),"h":ca.Sparsity.dense(2,2)},{"print_problem":True,"osqp.verbose":False})
    inputs = {"h": ca.DM([[1,0.2],[0.2,1]]),"g":ca.vertcat(1,2),"a":ca.horzcat(1,1),"lba":-1,"uba":1}
    extra_options = ["-Wno-endif-labels","-Wno-unused-variable","-Wno-strict-prototypes"]
    # No extra options on windows
    if os.name == 'nt':
      extra_options = []
    cg = self.check_codegen(c,inputs=c.convert_in(inputs),std="c99",extra_options=extra_options,extralibs=["osqp"],debug_mode=True)
    F = cg["F"].wrap()
    print("memful_external")
    self.check_codegen(F,inputs=c.convert_in(inputs),valgrind=True,main=True,extralibs=[cg["libname"]],debug_mode=True)

    # Check wrapper Function
    for X in [ca.MX,ca.SX]:

        g = X.sym("g",2)
        inputs["g"] = g

        f = ca.Function('wrapper',[g],[c(**inputs)['x']])
        print(f)
        # No extra options on windows
        if os.name == 'nt':
          extra_options = []
        cg = self.check_codegen(f,inputs=[ca.vertcat(1,2)],std="c99",extra_options=extra_options,extralibs=["osqp"],debug_mode=True)
        F = cg["F"].wrap()
        print("memful_external")
        self.check_codegen(F,inputs=[ca.vertcat(1,2)],valgrind=True,main=True,extralibs=[cg["libname"]],debug_mode=True)


  @requires_conic('osqp')
  def test_memful_jit(self):
    if not args.run_slow: return
    c = ca.conic("conic","osqp",{"a":ca.Sparsity.dense(1,2),"h":ca.Sparsity.dense(2,2)},{"print_problem":True,"osqp.verbose":False})
    inputs = {"h": ca.DM([[1,0.2],[0.2,1]]),"g":ca.vertcat(1,2),"a":ca.horzcat(1,1),"lba":-1,"uba":1}
    extra_options = ["-Wno-endif-labels","-Wno-unused-variable","-Wno-strict-prototypes"]
    # No extra options on windows
    if os.name == 'nt':
      extra_options = []

    # Check wrapper Function
    for X in [ca.MX,ca.SX]:
        g = X.sym("g",2)
        inputs["g"] = g

        libdir = ca.GlobalOptions.getCasadiPath()
        includedir = ca.GlobalOptions.getCasadiIncludePath()
        
        compiler_flags = extra_options
        linker_flags = []
        if os.name=='nt':
            compiler_flags.append("/I"+includedir)
            linker_flags.append("/LIBPATH:"+libdir)
            linker_flags.append("osqp.lib")
        else:
            compiler_flags.append("-I"+includedir)
            compiler_flags.append("-L"+libdir)
            compiler_flags.append("-losqp")
            linker_flags.append("-I"+includedir)
            linker_flags.append("-L"+libdir)
            linker_flags.append("-losqp")
            
        f = ca.Function('wrapper',[g],[c(**inputs)['x']])
        jit_opts = {"jit":True,"compiler":"shell","jit_options":{"compiler_flags":compiler_flags,"linker_flags": linker_flags,"verbose":True}}
        fj = ca.Function('wrapper',[g],[c(**inputs)['x']],jit_opts)
        self.checkfunction_light(f,fj,inputs=[ca.vertcat(1,2)])
        
        self.check_codegen(f,inputs=[ca.vertcat(1,2)],std="c99",external_opts=jit_opts,extra_options=extra_options,extralibs=["osqp"],debug_mode=True)
        
   
  def test_cse(self):
    for X in [ca.SX,ca.MX]:
        x = X.sym("x",2)
        y = X.sym("y",2)
        p = X.sym("p",2)
        z=p*ca.cos(x+5*y)
        z2=p*ca.cos(x+5*y)
        w = z-z2
        self.assertFalse(w.is_zero())
        res = ca.cse(w)
        self.assertTrue(res.is_zero())
        f1 = ca.Function('f',[x,y,p],[w],{"cse":False})
        f2 = ca.Function('f',[x,y,p],[w],{"cse":True})
        self.assertTrue(f1.n_instructions()>3)
        self.assertTrue(f2.n_instructions()<=3)

  def test_cse_call(self):
    for X in [ca.SX,ca.MX]:
        x = X.sym("x")
        f = ca.Function("f",[x],[x**2],["x"],["y"],{"never_inline":True})
        fcopy = ca.Function("f",[x],[x**2],["x"],["y"],{"never_inline":True})
 
        y = ca.cse(ca.vertcat(f(ca.sin(x)),fcopy(ca.sin(x)),f(ca.sin(x))))
        
        print(y)
        
        if X is ca.MX:
            self.assertTrue("vertcat(@1, @1, @1)" in str(y))
        else:
            self.assertTrue("[@1, @1, @1]" in str(y))

  @memory_heavy()
  def test_stop_diff(self):
    x = ca.MX.sym("x")
    y = ca.sin(x)
    
    max_order = 2
    for order in range(1, max_order+1):
    
        test_z = 3*x**(order-1)
        fa = ca.Function('f',[x],[ca.stop_diff(test_z,order)],['x'],['z'])
        fb = ca.Function('f',[x],[test_z],['x'],['z'])
        self.checkfunction(fa, fb, inputs=[ca.DM.rand(1,1)])
        
        z = ca.stop_diff(y,order)
        ca.DM.rng(1)
        
        v0 = ca.DM.rand(1)#v.sparsity())
        expr0 = ca.DM.rand(1)#expr.sparsity())
        
        for op in [ca.gradient,ca.jacobian,lambda expr,v : ca.jtimes(expr,v,v0),lambda expr,v : ca.jtimes(expr,v,expr0,True)]:
            yd = y
            zd = z
            for i in range(max_order+2):
                
                f = ca.Function('f',[x],[yd,zd])
                x0 = ca.DM.rand(f.sparsity_in(0))
                res = f(x0)
                

                if i>=order:
                    self.assertTrue(res[1].is_zero())
                else:
                    self.assertEqual(res[0],res[1])
        
                if i<order:
                    yd = op(yd,x)
                    zd = op(zd,x)
        
        f = ca.Function('f',[x],[z],['x'],['z'])
        f_normal = ca.Function('f',[x],[y],['x'],['z'])
        f_ref = ca.Function('f',[x],[test_z],['x'],['z'])
        
        for i in range(3):
            inputs = [ca.DM.rand(f.sparsity_in(i)) for i in range(f.n_in())]
            if i>=order:
                res = f.call(inputs,False,False)
                res_ref = f_ref.call(inputs,False,False)
                for e,e_ref in zip(res,res_ref):
                    if e_ref.is_zero():
                        self.assertTrue(e.is_zero())
                    else:
                        self.assertFalse(e.is_zero())
            else:
                self.checkfunction_light(f,f_normal,inputs=inputs)
                
            if i<order:
                
                in_labels = f.name_in()+['fwd:'+e for e in f.name_in()]+['adj:'+e for e in f.name_out()]
                out_labels = ['fwd:'+e for e in f.name_out()]+['adj:'+e for e in f.name_in()]
                for e_in in f.name_in():
                    for e_out in f.name_out():
                        out_labels.append('jac:%s:%s' % (e_out,e_in))
                        out_labels.append('grad:%s:%s' % (e_out,e_in))
                        
                f,f_normal,f_ref = [fe.factory('f',in_labels,out_labels) for fe in [f,f_normal,f_ref]]

  @memory_heavy()
  def test_stop_diff_cross(self):
    x = ca.MX.sym("x",2)
    y = ca.sin(x[0])+ca.exp(x[1])+ca.cos(x[0]**2*x[1]**2)
    
    max_order = 2
    for order in range(1, max_order+1):
    
        test_z = 3*x**(order-1)
        fa = ca.Function('f',[x],[ca.stop_diff(test_z,order)],['x'],['z'])
        fb = ca.Function('f',[x],[test_z],['x'],['z'])
        self.checkfunction(fa, fb, inputs=[ca.DM.rand(1,1)])
        
        z = ca.stop_diff(y,order)
        ca.DM.rng(1)

        for op in [ca.jacobian,lambda expr,v : ca.jtimes(expr,v,ca.DM.rand(v.sparsity())),lambda expr,v : ca.jtimes(expr,v,ca.DM.rand(expr.sparsity()),True)]:
            yd = y
            zd = z
            for i in range(max_order+1):
                
                f = ca.Function('f',[x],[yd,zd])
                x0 = ca.DM.rand(f.sparsity_in(0))
                res = f(x0)
                

                if i>=order:
                    self.assertTrue(res[1].is_zero())
                else:
                    self.checkarray(res[0],res[1])
        
                if i<order:
                    ca.DM.rng(5+i)
                    yd = op(yd,x)
                    ca.DM.rng(5+i)
                    zd = op(zd,x)
        
        f = ca.Function('f',[x],[z],['x'],['z'])
        f_normal = ca.Function('f',[x],[y],['x'],['z'])
        f_ref = ca.Function('f',[x],[test_z],['x'],['z'])
        
        if order<=1:
            for i in range(3):
                inputs = [ca.DM.rand(f.sparsity_in(i)) for i in range(f.n_in())]
                if i>=order:
                    res = f.call(inputs,False,False)
                    res_ref = f_ref.call(inputs,False,False)
                    for e,e_ref in zip(res,res_ref):
                        if e_ref.is_zero():
                            self.assertTrue(e.is_zero())
                        else:
                            self.assertFalse(e.is_zero())
                else:
                    self.checkfunction_light(f,f_normal,inputs=inputs)
                    
                if i<order:
                    
                    in_labels = f.name_in()+['fwd:'+e for e in f.name_in()]+['adj:'+e for e in f.name_out()]
                    out_labels = ['fwd:'+e for e in f.name_out()]+['adj:'+e for e in f.name_in()]
                    for e_in in f.name_in():
                        for e_out in f.name_out():
                            out_labels.append('jac:%s:%s' % (e_out,e_in))
                            
                    f,f_normal,f_ref = [fe.factory('f',in_labels,out_labels) for fe in [f,f_normal,f_ref]]

  @memory_heavy()
  def test_stop_diff_p(self):
    # Drop this -> not implemented
    x = ca.MX.sym("x")
    y = ca.MX.sym("y")
    yy = ca.sin(y*x)
    
    max_order = 2
    
    for order in range(1, max_order+1):
        ca.DM.rng(1)
        print("order",order)
        test_z = 3*(y*x)**(order-1)+y**order
        test_z = ca.sin(y)*x**(order-1)
        """fa = Function('f',[x,p],[stop_diff(test_z,x,order)],['x','p'],['z'])
        fb = Function('f',[x,p],[test_z],['x','p'],['z'])
        print()
        self.checkfunction(fa, fb, inputs=[DM.rand(1,1),DM.rand(1,1)])
        continue"""
        
        z = ca.stop_diff(yy,x,order)

        
        v0 = ca.DM.rand(1)#v.sparsity())
        expr0 = ca.DM.rand(1)#expr.sparsity())
        
        
        for op in [ca.gradient,ca.jacobian,lambda expr,v : ca.jtimes(expr,v,v0),lambda expr,v : ca.jtimes(expr,v,expr0,True)]:
            yd = yy
            zd = z
            for i in range(max_order+2):
                
                f = ca.Function('f',[x,y],[yd,zd])
                x0 = ca.DM.rand(f.sparsity_in(0))
                p0 = ca.DM.rand(f.sparsity_in(1))
                res = f(x0, p0)
                

                if i>=order:
                    self.assertTrue(res[1].is_zero())
                else:
                    self.assertEqual(res[0],res[1])
        
                if i<order:
                    yd = op(yd,x)
                    zd = op(zd,x)
        
        f = ca.Function('f',[x,y],[z],['x','y'],['z'])
        f_normal = ca.Function('f',[x,y],[yy],['x','y'],['z'])
        f_ref = ca.Function('f',[x,y],[test_z],['x','y'],['z'])
        
        for i in range(3):
            inputs = [ca.DM.rand(f.sparsity_in(i)) for i in range(f.n_in())]
            if i>=order:
                res = f.call(inputs,False,False)
                res_ref = f_ref.call(inputs,False,False)
                f_ref.disp(True)
                f.generate('f.c')
                for e,e_ref,name in zip(res,res_ref,f.name_out()):
                    print("check",name,e,e_ref)
                    if e_ref.is_zero():
                        self.assertTrue(e.is_zero())
                    else:
                        self.assertFalse(e.is_zero())
            else:
                self.checkfunction_light(f,f_normal,inputs=inputs)
                
            if i<order:
                
                in_labels = f.name_in()+['fwd:'+e for e in f.name_in()]*0+['adj:'+e for e in f.name_out()]
                out_labels = f.name_out()+['fwd:'+e for e in f.name_out()]*0+['adj:'+e for e in f.name_in()]
                #for e_in in f.name_in():
                #    for e_out in f.name_out():
                #        out_labels.append('jac:%s:%s' % (e_out,e_in))
                #        out_labels.append('grad:%s:%s' % (e_out,e_in))
                        
                f,f_normal,f_ref = [fe.factory('f',in_labels,out_labels) for fe in [f,f_normal,f_ref]]
            
        
  @requires_conic("qpoases")
  def test_no_hess2(self):
    y = ca.MX.sym("y",2)

    f = ca.Function("foo",[y],[ca.vertcat(ca.sin(y[0]*y[1]),ca.cos(y[0]*y[1]))],{"never_inline":True,"print_out":True})


    opti = ca.Opti()

    x = opti.variable(2)

    opti.subject_to(ca.no_hess(f(x**2))>=0)

    opti.minimize(ca.sumsqr(x-2))

    opti.solver("sqpmethod")

    F = opti.to_function("F",[x],[x])
    F.generate('F.c')
    with open('F.c','r') as inp:
        code = inp.read()
    for m in re.findall(r"\w+FS",code):
        self.assertTrue(len(m.split("_"))<=2)
    for m in re.findall(r"\w+foo",code):
        pass
        # bug 3019
        #self.assertTrue(len(m.split("_"))<=2)

  def test_issue_1522(self):
    V2 = ca.MX.sym("X",8)
    x,travel = ca.vertsplit(V2,[0,6,8])
    xs = ca.vertsplit(x,[0,3,6])
    travels = ca.vertsplit(travel,[0,1,2])

    dist = 0

    for j in range(2):
      dist+=ca.sum1((xs[0]-(xs[j]+travels[j]))**2)

    nlp = {"x":V2,"f":-dist}


    p = ca.MX(0,1)
    f = ca.Function('f',[V2,p],[-dist])


    for compact in [True,False]:
        for F in [f, f.expand()]:
            HF = F.reverse(1).jac_sparsity(0,0,compact)
            self.assertTrue(HF.is_symmetric())

  def test_issue_3074(self):
    x = ca.MX.sym("x",4)
    y = x[:2]
    y = 1
    f = ca.Function('f',[x],[ca.dot(y,y)])
    print(f)
    fr = f.reverse(1)
    print(fr)
    #print(fr.jac_sparsity())
  
    x = ca.MX.sym("x",ca.sparsify(ca.DM([0,1,1,0])).sparsity())
    y = ca.dot(x,x)
    f = ca.Function('f',[x],[ca.dot(y,y)])
    print(f)
    fr = f.reverse(1)
    print(fr)
    print(fr.jac_sparsity())
    
    data = "jhpnnagiieahaaaadaaaaaaaaaaaaaaaaafaegaakaaaaaaaneifgefhogdgehjgpgogbaaaaaaapaaaaaaaegfgmgbgjhpfbgchhgfhngfgogehdhaaaaaacaaaaaaahaaaaaaaaaaaaaaabababababababaaaaaaaaaaaaaaaaahaaaaaaaaaaaaaaaegfaaaaaaaaaaaaaaabaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaeggaaaaaaaaaaaaaaacaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaacaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabaaaaaaaaaaaaaaachbaaaaaaaaaaaaaaaegkaaaaaaaaaaaaaaagaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaagaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabaaaaaaaaaaaaaaacaaaaaaaaaaaaaaadaaaaaaaaaaaaaaaeaaaaaaaaaaaaaaafaaaaaaaaaaaaaaaegeaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaachdaaaaaaaaaaaaaaaegiaaaaaaaaaaaaaaaeaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaeaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabaaaaaaaaaaaaaaacaaaaaaaaaaaaaaadaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaahaaaaaaaaaaaaaaacaaaaaaajgadcaaaaaaajgbdcaaaaaaajgcdcaaaaaaajgddcaaaaaaajgedcaaaaaaajgfdcaaaaaaajggdaaaaaaaaaaaaaaaaaabagaaaaaaadhpgfhchdgfgbahaaaaaaakgjgehpfehngahaaaaaaaaaaaaaaaafaaaaaaadgmgbgoghgaaegbaaaaaaaaaaaaaaaaebababaaabababaaapbfilobfilobfnpdmfpicmfpicmfpnpdaaaaaeaaaaaaaaaaaaaaaabakdmiadcooijhfeodaaaaaaaaaaaaaaaabaaaaaaaocdaaaaaaangehihaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaachfaaaaaaaaaaaaaaahaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabaaaaaaahaaaaaaaaaaaaaaaegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaaeaaaaaaaehjgngfgegndaaaaaacaaaaaaaaaaaaaaaegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaagaaaaaaacgpgegjhocihegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaaiaaaaaaacgpgegjhocghpfihchbaaaaaaaaaaaaaaaegndaaaaaacaaaaaaaaaaaaaaaegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaalaaaaaaaegfgchiccgpgegjhocihjcegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaanaaaaaaaegfgchiccgpgegjhocghpfihjcchbaaaaaaaaaaaaaaaegndaaaaaagaaaaaaaaaaaaaaaegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaaiaaaaaaacgpgegjhocggpfihegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaaiaaaaaaacgpgegjhocbgpfihegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaajaaaaaaabgdgdgfgmgocbgpfihegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaakaaaaaaabgdgdgfgmgocngbgpfihegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaalaaaaaaabgdgdgfgmgoccgpfihocfhegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaalaaaaaaabgdgdgfgmgoccgpfihocjhchcaaaaaaaaaaaaaaaegmcaaaaaaadaaaaaaaaaaaaaaaachdaaaaaaaaaaaaaaaegmcaaaaaaadaaaaaaaaaaaaaaaachdaaaaaaaaaaaaaaaegndaaaaaaeaaaaaaaaaaaaaaaegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaagaaaaaaacgpgegjhochgegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaagaaaaaaacgpgegjhocdgegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaagaaaaaaacgpgegjhocngegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaalaaaaaaabgdgdgfgmgoccgpfihoccgcheaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaahaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaa"


    f = ca.StringDeserializer(data).unpack()
    print(f.reverse(1))
    f.reverse(1).jac_sparsity()

  def test_issue_4345(self):
    # SXFunction must match MXFunction adjoint sparsity; a structurally-zero
    # sensitivity (unused input) should come out 0nz, not full-pattern.
    for X in [ca.SX,ca.MX]:
      x = X.sym("x", ca.Sparsity.upper(3))
      y = X.sym("y", ca.Sparsity.upper(3))
      for z, adj_y_nnz in [(x**2,          0),   # y unused        -> adj_y 0nz
                           (x**2 + y[0,0], 6),   # y partially used -> adj_y 6nz
                           (x*y,           6)]:  # y fully used     -> adj_y 6nz
        adj = ca.Function("f", [x, y], [z]).reverse(1)
        self.assertEqual(adj.nnz_out(0), 6)            # adj_x: full input pattern
        self.assertEqual(adj.nnz_out(1), adj_y_nnz)    # adj_y: SX must agree with MX

  @requires_nlpsol("ipopt")
  def test_issue_3134(self):
    p = ca.MX.sym("p")
    v = ca.MX.sym("v")

    u = ca.MX.sym("u")
    x = ca.vertcat(p,v)

    rhs = ca.vertcat(v,u)

    t0 = ca.MX.sym("t0")
    DT = ca.MX.sym("DT")

    F = ca.Function('F', [x, u, t0, DT], [x + DT * rhs], ['x0', 'u', 't0', 'DT'], ['xf'])


    x = []
    g = []

    Xk = ca.MX.sym("Xk",2)
    x.append(Xk)

    T = ca.MX.sym("T")
    x.append(T)

    g.append(T)

    t0 = ca.MX.sym("t0")
    x.append(t0)

    # Note: (t0+T)-(t0+T/2) does not simplify; sparsity pattern indicates dependence on T and t0
    DTs = [T/2,(t0+T)-(t0+T/2)]

    for i in range(2):
        Uk = ca.MX.sym("Uk")
        x.append(Uk)
        Xk_next = ca.MX.sym("Xk_next",2)
        x.append(Xk_next)
        Xf = F(x0=Xk, u=Uk, t0=t0, DT=DTs[i])["xf"]
        g.append(Xk_next-Xf)
        Xk = Xk_next

    nlp = {"x": ca.vvcat(x),"f": T, "g": ca.vvcat(g)}

    solver = ca.nlpsol("solver","ipopt",nlp)
   
  def test_which_depends(self):
    x = ca.SX.sym("x",2)
    y = ca.SX.sym("y",2)
    
    z = x+y
    z2 = x*y
    z3 = ca.sin(ca.dot(x,y))
    z4 = ca.cos(x)
    z5 = ca.sin(y)
    z6 = 12
    z7 = ca.vertcat(x[0],y[1])
    z8 = ca.vertcat(x[0]**2,y[1]**2)
    
    f = ca.Function("f",[x,y],[z,z2,z3,z4,z5,z6,z7,z8],["x","y"],["z","z2","z3","z4","z5","z6","z7","z8"])
    
    def mychecks(f,exempt=False):
        self.assertEqual(f.which_depends("x",["z"]),[True,True])
        self.assertEqual(f.which_depends("x",["z"],2),[False,False])
        if not exempt:
            self.assertEqual(f.which_depends("x",["z2"],2),[False,False])
        self.assertEqual(f.which_depends("x",["z3"],2),[True,True])
        self.assertEqual(f.which_depends("x",["z4"]),[True,True])
        self.assertEqual(f.which_depends("x",["z4"],2),[True,True])
        self.assertEqual(f.which_depends("x",["z5"]),[False,False])
        self.assertEqual(f.which_depends("x",["z5"],2),[False,False])
        self.assertEqual(f.which_depends("x",["z5","z6"]),[False,False])
        if not exempt:
            self.assertEqual(f.which_depends("x",["z","z2","z3","z4","z5","z6"],1,True),[True, True]+[True, True]+ [True]+[True,True]+[False, False]+[False])
            self.assertEqual(f.which_depends("x",["z","z2","z3","z4","z5","z6"],2,True),[False, False]+[False, False]+ [True]+[True,True]+[False, False]+[False])
        self.assertEqual(f.which_depends("x",["z","z3","z4","z5","z6"],1,True),[True, True]+ [True]+[True,True]+[False, False]+[False])
        self.assertEqual(f.which_depends("x",["z","z3","z4","z5","z6"],2,True),[False, False]+ [True]+[True,True]+[False, False]+[False])
        self.assertEqual(f.which_depends("x",["z7"]),[True,False])
        self.assertEqual(f.which_depends("y",["z7"]),[False,True])
        self.assertEqual(f.which_depends("x",["z7"],2),[False,False])
        self.assertEqual(f.which_depends("y",["z7"],2),[False,False])
        self.assertEqual(f.which_depends("x",["z8"],2),[True,False])
        self.assertEqual(f.which_depends("y",["z8"],2),[False,True])
        self.assertEqual(f.which_depends("x",["z7","z8"],1,True),[True, False]+[True, False])
        self.assertEqual(f.which_depends("x",["z7","z8"],2,True),[False, False]+[True, False])
        self.assertEqual(f.which_depends("y",["z7","z8"],1,True),[False, True]+[False, True])
        self.assertEqual(f.which_depends("y",["z7","z8"],2,True),[False, False]+[False, True])
    
    mychecks(f)
    
    if not args.run_slow: return
    cg = self.check_codegen(f,inputs=[ca.DM.rand(f.sparsity_in(i)) for i in range(f.n_in())],with_jac_sparsity=True,with_forward=True)
    mychecks(cg["F"],exempt=True)
   
   
  def test_cat_input(self):

    x = ca.MX.sym("x",2)
    
    fmx = ca.Function('f',[ca.vertcat(ca.MX(2,1),x)],[ca.sin(x)],{"always_inline":True})
    print(fmx)
    
    X = ca.MX.sym("X",4)#Sparsity.lower(3))
    print(fmx(X))
    
    f_eval = ca.Function('f_eval',[X],[fmx(X)])
    print(f_eval([np.nan,np.nan,3,np.nan]))
    

  def test_reshape_input(self):

    x = ca.MX.sym("x",6,2)
    xsx = ca.SX.sym("x",6,2)
    
    fmx = ca.Function('f',[ca.reshape(x,3,4)],[ca.sin(x)],{"always_inline":True})
    fsx = ca.Function('f',[ca.reshape(xsx,3,4)],[ca.sin(xsx)])
    
    ca.DM.rng(1)
    inputs = [ca.DM.rand(3,4)]
    for input in inputs:
        self.checkfunction_light(fmx,fsx,inputs=[input])
    
 
 
  def test_sparsity_cast_input(self):
  
    X = ca.MX.sym("X",ca.Sparsity.lower(3))
    f = ca.Function("f",[X],[X**2],{"always_inline":True})
  
    sp = ca.Sparsity.lower(3)
    N = sp.nnz()
    
    x = ca.MX.sym("x",N)
    xsx = ca.SX.sym("x",N)
    
    fmx = ca.Function('f',[ca.sparsity_cast(x,sp)],[ca.sin(x)],{"always_inline":True})
    fsx = ca.Function('f',[ca.sparsity_cast(xsx,sp)],[ca.sin(xsx)])
    
    inputs = [ca.sparsify(ca.DM([[1,0,0],[3,7,0],[1,8,9]])), ca.sparsify(ca.DM([[1,0,0],[3,0,0],[0,8,9]])), ca.DM([[1,1.1,8],[3,1.2,3],[1,8,9]])]
    for input in inputs:
        self.checkfunction_light(fmx,fsx,inputs=[input])
    
    for X in [ca.MX.sym("X",ca.Sparsity.lower(3)), ca.MX.sym("X",ca.Sparsity.diag(3)), ca.MX.sym("X",3,3)]:
        fcmx = ca.Function('f',[X],[fmx(X)])
        fcsx = ca.Function('f',[X],[fsx(X)])
        
        self.checkfunction_light(fcmx,fcsx,inputs=[ca.DM.rand(X.sparsity())])

  def test_external(self):
    if not args.run_slow: return
    with self.assertInException("config failed"):
        self.compile_external("F","assets/externalfun1.c")
    
    for n_args in [3,7]:
        [externalfun1,libname] = self.compile_external("F","assets/externalfun1.c",{"config_args":["foo"]*n_args},debug_mode=True)
        print(libname)
        r = externalfun1(0)
        print(externalfun1.stats())
        self.assertEqual(n_args+1,int(r[0]))
        self.assertEqual(0,int(r[1]))

        x = ca.MX.sym("x")
        
        g = ca.Function('g',[x],[    externalfun1(ca.sin(x))[0]+externalfun1(ca.cos(x))[0] ])
        
        self.assertEqual(g(0.3),2*(n_args+1))
        self.check_codegen(g,inputs=[0.3],extralibs=[libname])
        
        self.check_serialize(externalfun1,inputs=[0.3])
        self.check_serialize(g,inputs=[0.3])
        
        externalfun1 = None
        g = None
        import gc
        gc.collect()
        
  def test_cache(self):
    x = ca.MX.sym("x")
    f = ca.Function('f',[x],[x**2])
    
    ff = f.forward(1)
    ff2 = f.forward(1)
    
    h = hash(ff)
    
    self.assertEqual(hash(ff),hash(ff2))
    
    ff = None
    ff2 = None
    
    import gc
    gc.collect()
    
    f2a = f.forward(2)
    f2b = f.forward(2)
    
    # Should be recreated because it was purged
    ff = f.forward(1)
    #self.assertNotEqual(hash(ff),h) # does not pass on Windows, but this test looks unreliable anyway
    self.assertEqual(hash(f2a),hash(f2b))
    
  def test_copy_elision(self):
    import casadi as ca
    
    backup = ca.GlobalOptions.getCopyElisionMinSize()
    
    for enabled in [False,True]:
    
      ca.GlobalOptions.setCopyElisionMinSize(2 if enabled else -1)

      x = ca.MX.sym("x",10)

      a = ca.MX.sym("a",5)

      bs = ca.MX.sym("b",10)

      b = (12*bs).monitor("b")
      y = 6*x


      [s0,s1] = ca.vertsplit(x,[0,7,10])
      [ss0,ss1] = ca.vertsplit(ca.sin(x),[0,5,10])

      z = y[2:7]+a+b[2:7]
      z = ca.sum1(z)
      z2 = ca.vertcat(x,a,b)*9

      f = ca.Function('f',[x,a,bs],[z,z2, x, ca.sqrt(s1), ss0+ss1])

      f.disp(True)

      f.generate('f.c')
      
      code = open('f.c','r').read()
      
      self.assertEqual("wr3" in code, enabled)
      
      lines = """
  casadi_int i;
  casadi_real *rr, w0, *w1=w+2, w2, *w3=w+8, *w4=w+18, *w5=w+28, *w6=w+33, *w7=w+38, *w8=w+43, *w9=w+68;
  const casadi_real *cs, *wr3, *wr4, *wr6, *wr9;
  /* #0: @0 = 0 */
  w0 = 0.;
  /* #1: @1 = ones(1x5) */
  casadi_fill(w1, 5, 1.);
  /* #2: @2 = 6 */
  w2 = 6.;
  /* #3: @3 = input[0][0] */
  wr3 = arg[0] ? arg[0] : casadi_zeros;
  /* #4: @4 = (@2*@3) */
  for (i=0, rr=w4, cs=wr3; i<10; ++i) (*rr++)  = (w2*(*cs++));
  /* #5: @5 = @4[2:7] */
  for (rr=w5, cs=w4+2; cs!=w4+7; cs+=1) *rr++ = *cs;
  /* #6: @6 = input[1][0] */
  wr6 = arg[1] ? arg[1] : casadi_zeros;
  /* #7: @5 = (@5+@6) */
  for (i=0, rr=w5, cs=wr6; i<5; ++i) (*rr++) += (*cs++);
  /* #8: @2 = 12 */
  w2 = 12.;
  /* #9: @4 = input[2][0] */
  wr4 = arg[2] ? arg[2] : casadi_zeros;
  /* #10: @4 = (@2*@4) */
  for (i=0, rr=w4, cs=wr4; i<10; ++i) (*rr++)  = (w2*(*cs++));
  /* #11: @4 = monitor(@4, b) */
  CASADI_PRINTF("b:\\n");
  casadi_print_canonical(casadi_s0, w4);
  CASADI_PRINTF("\\n");
  /* #12: @7 = @4[2:7] */
  for (rr=w7, cs=w4+2; cs!=w4+7; cs+=1) *rr++ = *cs;
  /* #13: @5 = (@5+@7) */
  for (i=0, rr=w5, cs=w7; i<5; ++i) (*rr++) += (*cs++);
  /* #14: @0 = mac(@1,@5,@0) */
  casadi_mtimes_dense(w1, 1, 5, w5, 1, (&w0), 0);
  /* #15: output[0][0] = @0 */
  if (res[0]) res[0][0] = w0;
  /* #16: @0 = 9 */
  w0 = 9.;
  /* #17: @8 = vertcat(@3, @6, @4) */
  rr=w8;
  for (i=0, cs=wr3; i<10; ++i) *rr++ = *cs++;
  for (i=0, cs=wr6; i<5; ++i) *rr++ = *cs++;
  for (i=0, cs=w4; i<10; ++i) *rr++ = *cs++;
  /* #18: @8 = (@0*@8) */
  for (i=0, rr=w8, cs=w8; i<25; ++i) (*rr++)  = (w0*(*cs++));
  /* #19: output[1][0] = @8 */
  casadi_copy(w8, 25, res[1]);
  /* #20: output[2][0] = @3 */
  casadi_copy(wr3, 10, res[2]);
  /* #21: {NULL, @9} = vertsplit(@3) */
  wr9 = wr3+7;
  /* #22: @9 = sqrt(@9) */
  for (i=0, rr=w9, cs=wr9; i<3; ++i) *rr++ = sqrt( *cs++ );
  /* #23: output[3][0] = @9 */
  casadi_copy(w9, 3, res[3]);
  /* #24: @3 = sin(@3) */
  for (i=0, rr=w3, cs=wr3; i<10; ++i) *rr++ = sin( *cs++ );
  /* #25: {@6, @1} = vertsplit(@3) */
  casadi_copy(w3, 5, w6);
  casadi_copy(w3+5, 5, w1);
  /* #26: @6 = (@6+@1) */
  for (i=0, rr=w6, cs=w1; i<5; ++i) (*rr++) += (*cs++);
  /* #27: output[4][0] = @6 */
  casadi_copy(w6, 5, res[4]);
  return 0;
"""
      if enabled:
        print(code)
        print(lines.strip())
        self.assertTrue(lines.strip() in code)
      self.check_codegen(f,inputs=[ca.DM.rand(f.sparsity_in(i)) for i in range(f.n_in())])
      
      x = ca.MX.sym("x",100,10)

      i = ca.MX.sym("i")
      y = x[4:90,i]
      f = ca.Function("f",[i,x],[ca.sin(y)+x[7,i]])

      f.disp(True)

      f.generate('f.c')
      
      code = open('f.c','r').read()
      
      self.assertEqual("wr" in code, enabled)
      
      lines = """
  casadi_int i;
  casadi_real *rr, w1, w2, w3, *w4=w+1003;
  const casadi_real *cs, *wr0, *wr4;
  /* #0: @0 = input[1][0] */
  wr0 = arg[1] ? arg[1] : casadi_zeros;
  /* #1: @1 = 100 */
  w1 = 100.;
  /* #2: @2 = input[0][0] */
  w2 = arg[0] ? arg[0][0] : 0;
  /* #3: @3 = floor(@2) */
  w3 = floor( w2 );
  /* #4: @1 = (@1*@3) */
  w1 *= w3;
  /* #5: @4 = @0[(4:90;@1)] */
  wr4 = wr0 + 4 + (int)w1;
  /* #6: @4 = sin(@4) */
  for (i=0, rr=w4, cs=wr4; i<86; ++i) *rr++ = sin( *cs++ );
  /* #7: @1 = 100 */
  w1 = 100.;
  /* #8: @2 = floor(@2) */
  w2 = floor( w2 );
  /* #9: @1 = (@1*@2) */
  w1 *= w2;
  /* #10: @2 = @0[(7;@1)] */
  i = 7 + (int) w1;
  w2 = i>=0 && i<1000 ? wr0[i] : casadi_nan;
  /* #11: @4 = (@4+@2) */
  for (i=0, rr=w4; i<86; ++i) (*rr++) += w2;
  /* #12: output[0][0] = @4 */
  casadi_copy(w4, 86, res[0]);
  return 0;
"""
      if enabled:
        print(code)
        print(lines.strip())
        self.assertTrue(lines.strip() in code)
      self.check_codegen(f,inputs=[3,ca.DM.rand(f.sparsity_in(1))],std="c99")
    ca.GlobalOptions.setCopyElisionMinSize(backup)
    
  @skip("simde" not in ca.CasadiMeta.feature_list())
  @memory_heavy()
  @requiresPlugin(ca.Importer,"shell")
  def test_blazing_spline(self):
    def test_points(knots_orig):
        import itertools
        selector = [lambda e: e[0]-0.1,lambda e: e[0],lambda e: (e[0]+e[1])/2,lambda e: e[-1],lambda e: e[-1]+0.1]
        for s in itertools.product(selector,repeat=len(knots_orig)):
            yield [e(k) for e,k in zip(s,knots_orig)]

    for N in [1,2,3,4,5]:
      for precompute_coeff in [False, True]:
        for precompute_grid in [False, True]:
          if precompute_coeff and N>3: continue
          for parametric_grid in [False, True]:
            if parametric_grid and precompute_coeff and N>3: continue

            def knots_expand(k):
                return [k[0]]*3 + k + [k[-1]]*3
  
            knots0 = [0,0.2,0.5,0.8,1]
            knots1 = [0,0.1,0.5,0.9,1]
            knots2 = [0,0.4,1]
            knots3 = [0,0.3,1]
            knots4 = [0,0.2,1]
            x0 = [0.4,0.5,0.6,0.7,0.8][:N]
  
            knots_orig = [knots0, knots1, knots2, knots3, knots4][:N]
            knots = [knots_expand(e) for e in knots_orig][:N]
  
  
            x=ca.MX.sym("x",N)
  
            nc = int(np.prod([len(k)-4 for k in knots]))
            ca.DM.rng(1)
            data = ca.DM.rand(nc)
            C = ca.MX.sym("C",nc,1)
            Y = ca.bspline(x,C,knots,[3]*N,1)
            F_ref = ca.Function('f',[x,C],[Y])
  
            if parametric_grid:
              knot_dims = [len(k) for k in knots]
              if os.name=='nt':
                  F = ca.blazing_spline("F",knot_dims,{"precompute_coeff": precompute_coeff, "precompute_grid": precompute_grid, "pedantic_mode_order": "ignore", "pedantic_mode_size": "ignore", "jit":True,"jit_options":{"flags": ["/I"+ca.GlobalOptions.getCasadiIncludePath()]}})
              else:
                  F = ca.blazing_spline("F",knot_dims,{"precompute_coeff": precompute_coeff, "precompute_grid": precompute_grid, "pedantic_mode_order": "ignore", "pedantic_mode_size": "ignore", "jit":True,"jit_options":{"flags": ["-I"+ca.GlobalOptions.getCasadiIncludePath(),"-g","-ffast-math","-march=native"]}})
  
              # Stacked knots vector for substitution
              knots_stacked = []
              for k in knots:
                  knots_stacked.extend(k)
              knots_data = ca.DM(knots_stacked)
              K = ca.MX.sym("K",len(knots_stacked))
  
              # Wrap: substitute fixed knots into parametric F to get (x,C)->f
              F_call = F(x, C, K)
              F_sub = ca.substitute([F_call[0]], [K], [knots_data])
              F_test = ca.Function('F_test',[x, C], F_sub)
            else:
              F_ref2 = ca.Function('f',[x,C],[Y,ca.jacobian(Y,x)])
              F_ref2([0]*N,data)
              print(knots)
  
              if os.name=='nt':
                  F = ca.blazing_spline("F",knots,{"precompute_coeff": precompute_coeff, "precompute_grid": precompute_grid, "pedantic_mode_order": "ignore", "pedantic_mode_size": "ignore", "jit":True,"jit_options":{"flags": ["/I"+ca.GlobalOptions.getCasadiIncludePath()]}})
              else:
                  F = ca.blazing_spline("F",knots,{"precompute_coeff": precompute_coeff, "precompute_grid": precompute_grid, "pedantic_mode_order": "ignore", "pedantic_mode_size": "ignore", "jit":True,"jit_options":{"flags": ["-I"+ca.GlobalOptions.getCasadiIncludePath(),"-g","-ffast-math","-march=native"]}})
              F_test = ca.Function('F_test',[x, C], [F(x, C)[0]])
  
            all_points = list(test_points(knots_orig))
  
            if N>3:
                random_base.seed(1)
                all_points = random_base.sample(all_points,125)
  
            for a in all_points:
                self.checkfunction_light(F_test,F_ref,inputs=[ca.vcat(a),data])
  
            if not parametric_grid:
              F2 = ca.Function.deserialize(F.serialize({"debug":True}))
              for a in all_points:
                  self.checkfunction_light(F,F2,inputs=[ca.vcat(a),data])
  
            self.check_codegen(F_test,inputs=[ca.vcat(a),data],std="c99")
  
            if parametric_grid:
              y_sin = ca.substitute(ca.sin(F(x,C,K)), K, knots_data)
            else:
              y_sin = ca.sin(F(x,C))
  
            f = ca.Function('f',[x,C],[y_sin,ca.jacobian(y_sin,x),ca.gradient(y_sin,x)])
            fcse = ca.Function('fcse',[x,C],[y_sin,ca.jacobian(y_sin,x),ca.gradient(y_sin,x)],{"cse":True})
            self.checkfunction_light(f,fcse,inputs=[ca.vcat(x0),data])
  
            f2 = ca.Function('f',[x,C],[y_sin,ca.jacobian(y_sin,x),ca.gradient(y_sin,x)]+list(ca.hessian(y_sin,x)))
            fcse2 = ca.Function('fcse',[x,C],[y_sin,ca.jacobian(y_sin,x),ca.gradient(y_sin,x)]+list(ca.hessian(y_sin,x)),{"cse":True})
            self.checkfunction_light(f2,fcse2,inputs=[ca.vcat(x0),data])
  
            with capture_stdout() as result:
                fcse2.disp(True)
  
            self.assertEqual(result[0].count(" F_der_der("),1)
            self.assertEqual(result[0].count(" F_der("),0)
            self.assertEqual(result[0].count(" F("),0)
  
            with capture_stdout() as result:
                fcse.disp(True)
            self.assertEqual(result[0].count(" F_der_der("),0)
            self.assertEqual(result[0].count(" F_der("),1)
            self.assertEqual(result[0].count(" F("),0)
  
            # extract_parametric: check extractable parameter-only subexpressions
            if parametric_grid:
              y_J = ca.jacobian(F(x,C,K)[0], x)
              _,_,parametric_exprs = ca.extract_parametric(y_J, ca.vertcat(C,K))
              if precompute_coeff or precompute_grid:
                # grid: inv depends only on K; coeff: dC depends on C and K
                self.assertTrue(len(parametric_exprs)>0)
                self.assertTrue(any(not e.is_zero() for e in parametric_exprs))
              else:
                # neither: no parameter-only subexpressions
                self.assertEqual(len(parametric_exprs), 0)
            else:
              y_J = ca.jacobian(F(x,C)[0], x)
              _,_,parametric_exprs = ca.extract_parametric(y_J, C)
              if precompute_coeff:
                # dC = derivative_coeff(C) is a C-only subexpression
                self.assertTrue(len(parametric_exprs)>0)
                self.assertTrue(any(not e.is_zero() for e in parametric_exprs))
              else:
                # non-parametric grid/none: no parameter-only subexpressions
                self.assertEqual(len(parametric_exprs), 0)
  
            # Jacobian comparison
            FJ_ref = F_ref.jacobian()
  
            if parametric_grid:
              FJ_param = F.jacobian()
              adj_f = ca.MX.sym("adj_f", 1, 1)
              fj_out = FJ_param(x, C, K, adj_f)
              fj_x = ca.substitute(fj_out[0], K, knots_data)
              fj_C = ca.substitute(fj_out[1], K, knots_data)
              FJ_test = ca.Function('FJ_test',[x, C, adj_f],[fj_x, fj_C])
            else:
              FJ = F.jacobian()
              FJ.generate("FJ.c",{"main":True})
              FJ.generate_in("FJ_in.txt",[ca.vcat(a),data,0])
              FJ_test = FJ
  
            for a in all_points:
                self.checkfunction_light(FJ_test,FJ_ref,inputs=[ca.vcat(a),data,0])
  
            if not parametric_grid:
              print("ref")
              F_ref2([0]*N,data)
  
              FJ = FJ.jacobian()
              FJ.generate("FH.c",{"main":True})
              FJ.generate_in("FH_in.txt",[ca.vcat(x0),data,0,0,0])
              print(FJ)
  
              FJ_ref = FJ_ref.jacobian()
              import time
              t0 = time.time()
              FJ_ref(ca.vcat(a),data,0,0,0)
              print("FJ_ref" ,time.time()-t0)
              t0 = time.time()
              FJ(ca.vcat(a),data,0,0,0)
              print("FJ" ,time.time()-t0)
              for a in all_points:
                  self.checkfunction_light(FJ,FJ_ref,inputs=[ca.vcat(a),data,0,0,0])
  
            if parametric_grid:
              # Test that changing knots gives different results
              knots_alt_orig = [[0,0.3,0.6,0.7,1],[0,0.2,0.4,0.8,1],[0,0.5,1],[0,0.4,1],[0,0.3,1]][:N]
              knots_alt = [knots_expand(e) for e in knots_alt_orig][:N]
              knots_alt_stacked = []
              for k in knots_alt:
                  knots_alt_stacked.extend(k)
              knots_alt_data = ca.DM(knots_alt_stacked)

              y_alt = ca.substitute(F(x, C, K)[0], K, knots_alt_data)
              F_alt = ca.Function('F_alt',[x, C],[y_alt])
              res1 = F_test(ca.vcat(x0), data)
              res2 = F_alt(ca.vcat(x0), data)
              self.assertFalse(np.allclose(float(res1), float(res2)),
                  "Parametric knots should produce different results with different knot vectors")

          # lookup_mode test: uniform knots so 'exact' is valid
          x_lu = ca.MX.sym("x",N)
          knots_uniform_orig = [
              list(np.linspace(0,1,5)),
              list(np.linspace(0,1,5)),
              list(np.linspace(0,1,3)),
              list(np.linspace(0,1,3)),
              list(np.linspace(0,1,3)),
          ][:N]
          knots_uniform = [[k[0]]*3 + k + [k[-1]]*3 for k in knots_uniform_orig]
          nc_u = int(np.prod([len(k)-4 for k in knots_uniform]))
          ca.DM.rng(1)
          data_u = ca.DM.rand(nc_u)
          C_u = ca.MX.sym("C",nc_u,1)
          Y_u = ca.bspline(x_lu,C_u,knots_uniform,[3]*N,1)
          F_ref_u = ca.Function('f',[x_lu,C_u],[Y_u])

          points_u = list(test_points(knots_uniform_orig))
          if N>3:
              random_base.seed(1)
              points_u = random_base.sample(points_u,125)

          for lookup_mode in ["linear","exact","binary"]:
            opts_lu = {"precompute_coeff": precompute_coeff,
                       "precompute_grid": precompute_grid,
                       "lookup_mode": [lookup_mode]*N,
                       "pedantic_mode_order": "ignore",
                       "pedantic_mode_size": "ignore",
                       "jit": True}
            if os.name=='nt':
                opts_lu["jit_options"] = {"flags": ["/I"+ca.GlobalOptions.getCasadiIncludePath()]}
            else:
                opts_lu["jit_options"] = {"flags": ["-I"+ca.GlobalOptions.getCasadiIncludePath(),"-g","-ffast-math","-march=native"]}
            F_lu = ca.blazing_spline("F",knots_uniform,opts_lu)
            F_lu_test = ca.Function('F_lu_test',[x_lu, C_u], [F_lu(x_lu, C_u)[0]])

            for a in points_u:
                self.checkfunction_light(F_lu_test,F_ref_u,inputs=[ca.vcat(a),data_u])

            F_lu2 = ca.Function.deserialize(F_lu.serialize({"debug":True}))
            for a in points_u:
                self.checkfunction_light(F_lu,F_lu2,inputs=[ca.vcat(a),data_u])

  @skip("simde" not in ca.CasadiMeta.feature_list())
  def test_blazing_spline_pedantic(self):
    # The pedantic checks fire in BlazingSplineFunction::init() based purely
    # on knot counts; they don't require evaluation or JIT. We pass plain
    # monotonic float lists; the spline math validity is irrelevant here.

    # n-4 = 8 (power of 2) -> pedantic_mode_size trips
    knots_size_bad_1d = [[float(i) for i in range(12)]]            # n=12 -> n-4=8

    # Order trip: counts [17, 13] are not non-decreasing
    # 17 -> n-4=13 (ok), 13 -> n-4=9 (ok), so size check stays silent.
    knots_order_bad_2d = [
        [float(i) for i in range(17)],
        [float(i) for i in range(13)],
    ]

    # Cumulative-product trip: counts [8, 8] -> extents [4, 4] (individually
    # below the >=8 threshold), but prefix product 4*4 = 16 (pow2) trips.
    knots_cumul_bad_2d = [
        [float(i) for i in range(8)],
        [float(i) for i in range(8)],
    ]

    # Clean: counts [13, 17] non-decreasing, n-4 in {9, 13}, neither pow2,
    # cumulative product 9*13 = 117 also not pow2.
    knots_clean_2d = [
        [float(i) for i in range(13)],
        [float(i) for i in range(17)],
    ]
    knots_clean_1d = [[float(i) for i in range(13)]]               # n-4=9

    # ----- assertInException: 'error' mode raises -----
    with self.assertInException("powers of 2"):
      ca.blazing_spline("F", knots_size_bad_1d,
                     {"pedantic_mode_size": "error"})

    with self.assertInException("not increasing"):
      ca.blazing_spline("F", knots_order_bad_2d,
                     {"pedantic_mode_order": "error",
                      "pedantic_mode_size":  "ignore"})

    # Message should pinpoint the zero-based dim index of the offender.
    with self.assertInException("dim 0 (zero-based)"):
      ca.blazing_spline("F", knots_size_bad_1d,
                     {"pedantic_mode_size": "error"})

    # Bad mode value errors out (only checked when there is an offender to
    # report; clean configs short-circuit before the mode is interpreted).
    with self.assertInException("'pedantic_mode_size' must be one of"):
      ca.blazing_spline("F", knots_size_bad_1d,
                     {"pedantic_mode_size": "bogus"})

    # Cumulative-product trip: individuals are below the >=8 threshold but
    # the prefix product hits a power of 2.
    with self.assertInException("prefix product"):
      ca.blazing_spline("F", knots_cumul_bad_2d,
                     {"pedantic_mode_size": "error"})

    # ----- capture_stdout: 'warn' mode emits a warning, no raise -----
    with capture_stdout() as result:
      ca.blazing_spline("F", knots_size_bad_1d,
                     {"pedantic_mode_size": "warn"})
    self.assertTrue("powers of 2" in result[0] or "powers of 2" in result[1])

    with capture_stdout() as result:                               # order default = 'warn'
      ca.blazing_spline("F", knots_order_bad_2d,
                     {"pedantic_mode_size": "ignore"})
    self.assertTrue("not increasing" in result[0] or "not increasing" in result[1])

    # 'ignore' produces no warning text and no exception.
    with capture_stdout() as result:
      ca.blazing_spline("F", knots_size_bad_1d,
                     {"pedantic_mode_size":  "ignore",
                      "pedantic_mode_order": "ignore"})
    self.assertFalse("powers of 2" in result[0] or "powers of 2" in result[1])

    # ----- Clean configs must not raise even with 'error' on both knobs -----
    ca.blazing_spline("F", knots_clean_1d,
                   {"pedantic_mode_order": "error",
                    "pedantic_mode_size":  "error"})
    ca.blazing_spline("F", knots_clean_2d,
                   {"pedantic_mode_order": "error",
                    "pedantic_mode_size":  "error"})
    # Single dimension trivially satisfies the order constraint.
    ca.blazing_spline("F", knots_clean_1d, {"pedantic_mode_order": "error"})

    # ----- (n_knots - 5) path: precompute_coeff_=True && diff_order_>=1.
    # Only the jacobian child sees this, so the parent stays silent and
    # the trip fires when F.jacobian() builds the child. pedantic_mode_*
    # propagates from parent to child automatically; no jacobian_options
    # plumbing needed.
    # 13 knots -> n-4=9 (ok), n-5=8 (pow2), n-6=7 (ok).
    knots_n5_bad_1d = [[float(i) for i in range(13)]]

    # Default pedantic_mode_size = 'error' propagates -> jacobian raises.
    with self.assertInException("(n_knots - 5)"):
      F = ca.blazing_spline("F", knots_n5_bad_1d)
      F.jacobian()
    # Message also reports the diff order.
    with self.assertInException("diff order 1"):
      F = ca.blazing_spline("F", knots_n5_bad_1d)
      F.jacobian()

    # Parent 'warn' propagates -> child warns, no exception.
    with capture_stdout() as result:
      F = ca.blazing_spline("F", knots_n5_bad_1d,
                         {"pedantic_mode_size": "warn"})
      F.jacobian()
    self.assertTrue("(n_knots - 5)" in result[0] or "(n_knots - 5)" in result[1])

    # Parent 'ignore' propagates -> child silent too (no exception, no log).
    with capture_stdout() as result:
      F = ca.blazing_spline("F", knots_n5_bad_1d,
                         {"pedantic_mode_size": "ignore"})
      F.jacobian()
    self.assertFalse("(n_knots - 5)" in result[0] or "(n_knots - 5)" in result[1])

    # jacobian_options still wins when explicitly provided.
    with capture_stdout() as result:
      F = ca.blazing_spline("F", knots_n5_bad_1d,
                         {"jacobian_options": {"pedantic_mode_size": "ignore"}})
      F.jacobian()
    self.assertFalse("(n_knots - 5)" in result[0] or "(n_knots - 5)" in result[1])

    # 15 knots -> n-4=11, n-5=10, n-6=9 all clean at every diff order.
    knots_jac_clean_1d = [[float(i) for i in range(15)]]
    F = ca.blazing_spline("F", knots_jac_clean_1d,
                       {"pedantic_mode_size": "error"})
    F.jacobian()

    # ----- (n_knots - 6) path: !precompute_coeff_ && diff_order_>=2.
    # Reached via a grandchild built by F.jacobian().jacobian().
    # 14 knots -> n-4=10 (ok), n-5=9 (ok), n-6=8 (pow2).
    knots_n6_bad_1d = [[float(i) for i in range(14)]]

    # Default 'error' propagates two levels deep.
    with self.assertInException("(n_knots - 6)"):
      F = ca.blazing_spline("F", knots_n6_bad_1d,
                         {"precompute_coeff": False})
      F.jacobian().jacobian()
    with self.assertInException("diff order 2"):
      F = ca.blazing_spline("F", knots_n6_bad_1d,
                         {"precompute_coeff": False})
      F.jacobian().jacobian()

  def test_noncanonical_sparsity(self):
    x = ca.MX.sym("x",4,4)
    y = ca.MX.sym("y")
    
    f = ca.Function("f",[x,y],[3*8])
    
    self.check_codegen(f,inputs=[ca.DM.rand(4,4),1],opts={"force_canonical":False})
    self.check_codegen(f,inputs=[ca.DM.rand(4,4),1],opts={"force_canonical":True})
    
  def test_options_sanitize(self):
      def canonical(e):
        if isinstance(e,dict):
            return [(k,canonical(e[k])) for k in sorted(e.keys())]
        return e
      def cmp(a,b):
        assert str(canonical(a))==str(canonical(b))
      res = ca.Options.sanitize({"foo": {"bar": {"baz.goot": 7}}})
      ref = {"foo":{"bar": {"baz": {"goot":7}}}}
      cmp(res,ref)
      res = ca.Options.sanitize({"foo.baz": 1, "foo":{"bar": 9}})
      ref = {'foo': {'bar': 9, 'baz': 1}}
      cmp(res,ref)
      res = ca.Options.sanitize({"foo.baz.bar": 1, "foo":{"baz": {"go":8}}})
      ref = {"foo":{"baz":{"bar": 1, "go": 8}}}
      cmp(res,ref)
      res = ca.Options.sanitize({"foo":{"bar": 9}, "foo.baz": 1})
      ref = {"foo":{"bar":9, "baz":1}}
      cmp(res,ref)
      res = ca.Options.sanitize({"foo.bar.baz": 1, "foo.bar.mat":3,"foo.lib":5, "foo.bar.qint.pad":5,"foo":{"bar": {"got": 9}},"bar":{"abc.def":7}})
      ref = {'bar': {'abc': {'def': 7}}, 'foo': {'bar': {'baz': 1, 'mat': 3, 'got': 9, 'qint': {'pad': 5}}, 'lib': 5}}
      cmp(res,ref)
      
      with self.assertInException("update_dict error"):
          res = ca.Options.sanitize({"foo": 9, "foo.baz": 1})
      with self.assertInException("update_dict error"):
          res = ca.Options.sanitize({"foo": {"baz":{"w": 5}}, "foo.baz": 1})
          
      res = ca.Options.sanitize({"foo": 1, "bar": None})
      print(res)
      ref = {"foo": 1}
      cmp(res,ref)
      # It is important that None/null entries are preserved deeper than toplevel
      # Use for augmented_options machinery in integrator
      res = ca.Options.sanitize({"foo": 1, "bar": {"r": None}})
      print(res)
      ref = {"foo": 1, "bar": {"r": None}}
      cmp(res,ref)
      
  def test_simplify(self):
    x = ca.MX.sym("x",5)
    y = ca.MX.sym("y",5)
    f = ca.Function('test_func', [x, y], [x + y])

    # No name argument: default simplify flow, original name preserved
    f_simplified = f.transform()
    self.assertIsInstance(f_simplified, ca.Function)
    self.assertEqual(f_simplified.name(), 'test_func')

    # Same name argument
    f_same_name = f.transform("test_func")
    self.assertIsInstance(f_same_name, ca.Function)
    self.assertEqual(f_same_name.name(), 'test_func')

    # Different name argument renames the result
    f_renamed = f.transform("new_name")
    self.assertIsInstance(f_renamed, ca.Function)
    self.assertEqual(f_renamed.name(), 'new_name')

    # Test that simplified functions work correctly
    self.checkfunction_light(f,f_simplified,inputs=[2.0,3.0])

    x = ca.MX.sym("x",5)
    y = ca.MX.sym("y",5)
    f = ca.Function('test_func', [x, y], [x])
    f_simplified = f.transform([["simplify","empty_inputs"]])

    self.assertEqual(f.nnz_in(0), 5)
    self.assertEqual(f_simplified.nnz_in(0), 5)   # x is used
    self.assertEqual(f_simplified.nnz_in(1), 0)   # y is unused -> emptied

    # Test that simplified functions work correctly
    self.checkfunction_light(f,f_simplified,inputs=[2.0,3.0])




  def test_simplify_ref_count(self):
  
    # getnonzeros, vertsplit

  
    def tests():
    
        for X in [ca.SX,ca.MX]:

            a = X.sym("a")
            b = X.sym("b")
            
            z = -(a-b) # Will get simplified
            
            yield ca.Function("f",[a,b],[z]), ca.Function("f",[a,b],[b-a])
            
            e = (a-b)
            z = -e     # Will not get simplified
            
            yield ca.Function("f",[a,b],[3*e,-e]), ca.Function("f",[a,b],[3*e,z])
            
                    
        x = ca.MX.sym("x",ca.Sparsity.lower(3))

        y = ca.project(x, ca.Sparsity.upper(3))

        yield ca.Function("f",[x],[y.nz[5]]), ca.Function("f",[x],[x.nz[5]])

        x = ca.MX.sym("x",ca.Sparsity.lower(3))

        y = ca.project(x, ca.Sparsity.upper(3))

        yield ca.Function("f",[x],[y[:,-1]]), ca.Function("f",[x],[x.nzref(ca.Sparsity.dense(3,1),[-1,-1,5])])
        
        return


        x = ca.MX.sym("x")
        
        z = ca.vertcat(x,ca.sin(x),ca.cos(x))
        
        y = z[:2]
        
        yield ca.Function("f",[x],[y]), ca.Function("f",[x],[ca.vertcat(x,ca.sin(x))])


    for f, ref in tests():
        ref.disp(True)
        fs = f.transform()
        
        fs.disp(True)
        
        fs.generate("f.c")
        ref.generate("ref.c")
        self.checkfunction_identical(fs,ref)
        
        ca.DM.rng(1)
        inputs = [ca.DM.rand(f.sparsity_in(i)) for i in range(f.n_in())]
        self.checkfunction_light(f,fs,inputs=inputs)
    
  def test_simplify_const_folding(self):

    ca.DM.rng(1)


    A = ca.MX(ca.DM.rand(2,2))
    B = ca.MX(ca.DM.rand(2,2))
    C = A @ B

    C = C @ C

    x = ca.MX.sym("x")
    f = ca.Function('f',[x],[(x*(2*C)+(2*C))*C[0]])


    fs = f.transform()
    inputs = [ca.DM.rand(f.sparsity_in(i)) for i in range(f.n_in())]
    self.checkfunction_light(f,fs,inputs=inputs)
    
    fs.disp(True)
    self.assertEqual(fs.n_nodes(),7)

  def test_simplify_combine_terms(self):
    # issue #4227: expression-level simplify routes through Function.simplify
    x = ca.SX.sym("x")
    y = ca.SX.sym("y")

    def numeq(expr, ref):
      fe = ca.Function("fe", [x, y], [expr])
      fr = ca.Function("fr", [x, y], [ref])
      for vx, vy in [(0.3, 1.7), (-2.0, 4.5)]:
        self.checkarray(fe(vx, vy), fr(vx, vy))

    def nnodes(expr):
      return ca.Function("f", [x, y], [expr]).n_nodes()

    # like terms get combined into a single weighted term
    s1 = ca.simplify(8*y-3*y)
    numeq(s1, 5*y)
    self.assertTrue(nnodes(s1) < nnodes(8*y-3*y))

    s2 = ca.simplify(8*(5*x)-3*(x+3+((12*x)-x)-17*x))
    numeq(s2, 55*x-9)
    self.assertTrue(nnodes(s2) < nnodes(8*(5*x)-3*(x+3+((12*x)-x)-17*x)))

    # the same combining is available as a transform pass on a Function
    f = ca.Function("f",[x],[8*x-3*x])
    fc = f.transform([["simplify","combine_terms"]])
    self.checkfunction_light(f,fc,inputs=[2.0])
    self.assertTrue(fc.n_nodes() < f.n_nodes())

    # free simplify(MX) is a no-op (unchanged, long-standing behaviour)
    a = ca.MX.sym("a")
    b = ca.MX.sym("b")
    e = -(a-b)
    self.assertTrue(ca.is_equal(ca.simplify(e), e, 2))

  def test_transform(self):
    # issue #4227: ordered transformation passes
    x = ca.SX.sym("x")
    y = ca.SX.sym("y")

    # free-function transform on an expression, "simplify" verb + tasks
    s = ca.transform(8*y-3*y, [["simplify","combine_terms"]])
    fe = ca.Function("fe",[x,y],[s])
    fr = ca.Function("fr",[x,y],[5*y])
    for vx,vy in [(0.3,1.7),(-2.0,4.5)]:
      self.checkarray(fe(vx,vy), fr(vx,vy))

    # free-function transform on a vector of expressions, transformed jointly
    sv = ca.transform([8*y-3*y, x+x], [["simplify","combine_terms"]])
    self.assertEqual(len(sv), 2)
    fev = ca.Function("fev",[x,y],sv)
    frv = ca.Function("frv",[x,y],[5*y, 2*x])
    for vx,vy in [(0.3,1.7),(-2.0,4.5)]:
      self.checkarray(fev(vx,vy)[0], frv(vx,vy)[0])
      self.checkarray(fev(vx,vy)[1], frv(vx,vy)[1])
    # the dict form is also available for the vector signature
    self.assertEqual(len(ca.transform([8*y-3*y, x+x], dict(combine_terms=True))), 2)

    # Function.transform with an ordered list; 0 = run a pass until fixed point
    a = ca.MX.sym("a")
    b = ca.MX.sym("b")
    f = ca.Function("f",[a,b],[-(a-b)])
    ft = f.transform([["simplify",0,"ref_count","const_folding"]])
    self.checkfunction_identical(ft, ca.Function("f",[a,b],[b-a]))
    self.assertTrue(ft.n_nodes() < f.n_nodes())

    # a positive integer runs a pass a fixed number of times
    ft2 = f.transform([["simplify",2,"ref_count"]])
    self.checkfunction_identical(ft2, ca.Function("f",[a,b],[b-a]))

    # name forms: default keeps the original name, explicit name renames
    self.assertEqual(f.transform().name(), "f")
    self.assertEqual(f.transform([["simplify","cse"]]).name(), "f")
    self.assertEqual(f.transform("g").name(), "g")
    self.assertEqual(f.transform("g", [["simplify","cse"]]).name(), "g")
    self.assertEqual(f.transform("g", dict(cse=True)).name(), "g")
    fr = f.transform("g", [["simplify",0,"ref_count","const_folding"]])
    self.checkfunction_identical(fr, ca.Function("g",[a,b],[b-a]))

    # no-argument, empty dict, and all-default booleans apply the default flow
    fd_noarg = f.transform()
    fd_empty = f.transform({})
    fd_bool = f.transform(dict(combine_terms=False))
    for fd in [fd_noarg, fd_empty, fd_bool]:
      self.checkfunction_identical(fd, ca.Function("f",[a,b],[b-a]))

    # an explicit empty pass list is a no-op (distinct from default flow)
    fd_none = f.transform([])
    self.checkfunction_identical(fd_none, f)

    # boolean shorthands forward to a single "simplify" pass
    s = ca.SX.sym("s")
    g = ca.Function("g",[s],[8*s-3*s])
    gc = g.transform(dict(combine_terms=True))
    self.checkfunction_light(g, gc, inputs=[2.0])
    self.assertTrue(gc.n_nodes() < g.n_nodes())
    with self.assertRaises(Exception):
      g.transform(dict(bogus=True))
    # the (passes, opts) form rejects boolean simplify options in opts
    with self.assertRaises(Exception):
      g.transform([["simplify","cse"]], dict(cse=True))

    # expand pass turns an MXFunction into an SXFunction
    fx = f.transform([["expand"]])
    self.assertEqual(fx.class_name(), "SXFunction")
    self.checkfunction_light(f, fx, inputs=[1.0, 2.0])

    # unknown verb / task raise
    with self.assertRaises(Exception):
      f.transform([["bogus"]])
    with self.assertRaises(Exception):
      f.transform([["simplify","nope"]])
    # a negative run count is rejected
    with self.assertRaises(Exception):
      f.transform([["simplify",-1,"ref_count"]])

    # external pass adopts external_transform, with an optional opts dict
    if sys.platform != 'darwin':
      libcasadi = ca.CasadiMeta.shared_library_prefix()+"casadi"
      g = ca.Function("g",[x],[x**2])
      with self.assertOutputs(["Doing","bar"],["Warning"]):
        gt = g.transform([["external", libcasadi,
                           "external_transform_test_success", {"foo": "bar"}]])
      self.checkfunction(g, gt, inputs=[3])

    # verbose prints the copy-pasteable passes and per-pass stats
    with self.assertOutputs(["transform(", "n_nodes"],[]):
      f.transform([["simplify",0,"ref_count"]], dict(verbose=True))
    with self.assertOutputs(["transform(", "n_nodes"],[]):
      f.transform(dict(verbose=True))

  def test_duplicate_check(self):
    x = ca.MX()
    f = ca.Function('f',[x],[2*x])
    self.assertEqual(f.numel_in(0), 0)
    f = ca.Function('f',[ca.MX.zeros(0,0)],[2*x])
    self.assertEqual(f.numel_in(0), 0)
    f = ca.Function('f',[ca.MX.ones(0,0)],[2*x])
    self.assertEqual(f.numel_in(0), 0)
    p = ca.MX.sym("p")
    with self.assertInException("input arguments must be purely symbolic"):
        f = ca.Function('f',[x.nz[p]],[2*x])
        self.assertEqual(f.numel_in(0), 0)
    x.nz[p] = 3.7
    with self.assertInException("input arguments must be purely symbolic"):
        f = ca.Function('f',[x],[3*x])
        self.assertEqual(f.numel_in(0), 0)

  def test_finite_diff(self):
    x = ca.MX.sym("x")
    for method in ["forward","backward","central","smoothing"]:
        f = ca.Function("f",[x],[ca.sin(x)],{"enable_fd":True,"enable_forward":False,"enable_reverse":False,"fd_method":method})
        fref = ca.Function("f",[x],[ca.sin(x)])
        J = f.jacobian()
        Jref = fref.jacobian()

        inputs = [1.1,0]
        self.checkfunction_light(J,Jref,inputs=inputs,digits=5)
        
        if args.run_slow:
          digits = 15
          if method=="smoothing" and os.name == 'nt':
              digits = 10# h is determined with pow(abstol,1/3) 
          self.check_codegen(J,inputs=inputs,std="c99",digits=digits)
      
       
  def test_is_diff_fd(self):

    ca.DM.set_precision(16)
    
    
    for n in [1,2]:
        for m in [1,3]:
            x = ca.MX.sym("x",n,m)
            y = ca.MX.sym("y",n,m)

            expr = [2*ca.sumsqr(x)*ca.sumsqr(y),y**2*ca.sumsqr(x)]
            fref = ca.Function("f",[x,y],expr,["x","y"],["a","b"],{"is_diff_in":[True,False]})
            f = ca.Function("f",[x,y],expr,["x","y"],["a","b"],{"is_diff_in":[True,False],"enable_fd":True,"print_in":True,"print_out":False,"enable_forward":False,"enable_reverse":False,"fd_method":"forward"})
    
            Jref = fref.jacobian()
            J = f.jacobian()
            self.assertTrue(Jref.is_diff_in(0))
            self.assertFalse(Jref.is_diff_in(1))

            self.assertEqual(Jref.nnz_in(2),0)
            self.assertEqual(J.nnz_in(2),0)
            
            Ffref = fref.forward(2)
            Ff = f.forward(2)
            
            self.assertEqual(Ffref.nnz_in(2),0)
            self.assertEqual(Ff.nnz_in(2),1)
            
            for F,Fref in [(J,Jref),(Ff,Ffref)]:
            
                print("F = ", F)
                print("Fref = ", Fref)

                
                self.assertEqual(F.n_in(), Fref.n_in())
                self.assertEqual(F.n_out(), Fref.n_out())
                 
                for i in range(F.n_in()):
                    self.assertEqual(F.name_in(i),Fref.name_in(i))
                    if not F.name_in(i).startswith("out_"):
                        self.assertEqual(F.sparsity_in(i),Fref.sparsity_in(i))
                    self.assertEqual(F.is_diff_in(i),Fref.is_diff_in(i))

                ca.DM.rng(1)
                inputs = [ca.DM.rand(F.sparsity_in(i)) for i in range(f.n_in())]
                inputs += f.call(inputs)
                if "fwd" in F.name():
                    inputs += [ca.DM.rand(F.sparsity_in(i)) for i in range(f.n_in())]
                if "adj" in F.name():
                    inputs += [ca.DM.rand(F.sparsity_out(i)) for i in range(f.n_out())]
                self.checkfunction_light(F,Fref,inputs=inputs,digits=5)
                
                if args.run_slow:
                    digits = 15
                    if sys.platform=="darwin":
                        # very sensitive to changes in fused-multiply-add versus non-fused
                        # Can be steered with -ffp-contract, but dind't find a universally working incantation
                        digits = 5
                    self.check_codegen(F,inputs=inputs,std="c99",digits=digits)
              
                inputs = [2,3] + [0,0]
                if "fwd" in F.name():
                    inputs += [ca.DM.rand(F.sparsity_in(i)) for i in range(f.n_in())]
                if "adj" in F.name():
                    inputs += [ca.DM.rand(F.sparsity_out(i)) for i in range(f.n_out())]
                    
                with capture_stdout() as result:
                    F(*inputs)
                self.assertTrue("2.00" in result[0])
                self.assertFalse("3.00" in result[0])
                
    ca.DM.set_precision(6)

  def test_activity(self):
    # issue #3019: value-level "signal activity" propagation (Function::activity).
    # True/active = possibly nonzero. The expected output activity is never hard-coded:
    # it is derived empirically by evaluating the SAME function on DMs whose active
    # input nonzeros are set to (distinct) nonzeros and inactive ones to 0, and taking
    # the union of nonzero outputs over a few value assignments (so the comparison is
    # both sound -- never claim zero when it can be nonzero -- and tight).
    def emp_activity(f, mask):
        emp = [False]*f.nnz_out()
        for t in range(1, 6):
            ins = []; off = 0
            for i in range(f.n_in()):
                ni = f.nnz_in(i)
                vals = [float(t*(off+j+1)) if mask[off+j] else 0.0 for j in range(ni)]
                ins.append(ca.DM(f.sparsity_in(i), vals)); off += ni
            flat = []
            for o in f.call(ins): flat += list(o.nonzeros())
            emp = [emp[k] or (flat[k] != 0) for k in range(f.nnz_out())]
        return emp
    def all_masks(n):
        return [[bool(b & (1 << i)) for i in range(n)] for b in range(2**n)]
    def check(f, masks, minimal=True):
        # activity must be SOUND (never predict zero where the value can be nonzero);
        # with minimal=True it must additionally be tight (exactly equal to empirical).
        for m in masks:
            emp = emp_activity(f, m); got = list(f.activity(m))
            for k in range(f.nnz_out()):
                self.assertTrue(got[k] or not emp[k],
                                "UNSOUND: %s mask %s output %d: activity says inactive but "
                                "empirically nonzero" % (f.name(), m, k))
            if minimal:
                self.assertEqual(got, emp, "not minimal: %s mask %s" % (f.name(), m))

    # Elementary ops, MX and SX, exhaustive over the input activity space
    for sym in [ca.MX.sym, ca.SX.sym]:
        x = sym("x"); y = sym("y")
        check(ca.Function("mul", [x, y], [x*y]),   all_masks(2))   # 0-annihilation
        check(ca.Function("add", [x, y], [x+y]),   all_masks(2))
        check(ca.Function("cs",  [x],    [ca.cos(x)]), all_masks(1))  # cos(0)=1: active from inactive
        check(ca.Function("sn",  [x],    [ca.sin(x)]), all_masks(1))  # sin(0)=0
        check(ca.Function("z",   [x],    [0*x]),    all_masks(1))  # constant zero

    # Partially-zero constant: per-element precision, not a whole-matrix is_zero() collapse
    x2 = ca.MX.sym("x2", 2)
    check(ca.Function("pc", [x2], [x2*ca.densify(ca.MX(ca.DM([0, 5])))]), all_masks(2))

    # mtimes: minimally conservative (a term is active only iff BOTH factors are); a
    # dependency analysis combining factors with OR would over-report e.g. [T,T,F,F].
    xr = ca.MX.sym("xr", 1, 2); yc = ca.MX.sym("yc", 2, 1)
    check(ca.Function("mt", [xr, yc], [ca.mtimes(xr, yc)]), all_masks(4))

    # solve: a zero rhs gives a zero result (X = A\0 = 0). A dense A couples the whole
    # column, so activity is tight (minimal); a structured (e.g. diagonal) A does not
    # couple, so the whole-column rule is only sound, not minimal -- but must never lie.
    A = ca.MX.sym("A", 2, 2); b = ca.MX.sym("b", 2, 1)
    check(ca.Function("slv", [A, b], [ca.solve(A, b)]), [[True]*4 + m for m in all_masks(2)])
    ad = ca.MX.sym("ad", 2); bd = ca.MX.sym("bd", 2, 1)
    check(ca.Function("sd", [ad, bd], [ca.solve(ca.diag(ad), bd)]),
          [[True]*2 + m for m in all_masks(2)], minimal=False)

    # The motivating case: 2nd-order sensitivity, exhaustive 256-mask sweep. In
    # particular the mask with zero forward seeds yields an all-inactive output.
    y = ca.MX.sym("y", 2)
    foo = ca.Function("foo", [y], [ca.vertcat(ca.sin(y[0]*y[1]), ca.cos(y[0]*y[1]))])
    fadj = foo.reverse(1).forward(1)
    check(fadj, all_masks(fadj.nnz_in()))

    # End-to-end: the inactive 2nd-order call node is dropped at construction (Call::create)
    args = []
    for i in range(fadj.n_in()):
        if fadj.name_in(i).startswith("fwd_"):
            args.append(ca.MX(fadj.size1_in(i), fadj.size2_in(i)))    # zero seed
        else:
            args.append(ca.MX.sym("a%d" % i, fadj.sparsity_in(i)))
    self.assertTrue(fadj.call(args)[0].is_zero())                  # removed -> structural zero
    args_nz = [ca.MX.sym("b%d" % i, fadj.sparsity_in(i)) for i in range(fadj.n_in())]
    self.assertFalse(fadj.call(args_nz)[0].is_zero())              # nonzero seeds -> call kept

    # no_hess Hessian-of-Lagrangian no longer embeds the 2nd-order derivative call
    x = ca.MX.sym("x", 2)
    lam = ca.MX.sym("lam", 2); sig = ca.MX.sym("sig")
    L = sig*ca.sumsqr(x-2) + ca.dot(lam, ca.no_hess(foo(x**2)))
    H = ca.hessian(L, x)[0]
    Hf = ca.Function("H", [x, lam, sig], [H])
    self.assertEqual(len(Hf.get_function()), 0)                   # no embedded call nodes remain
    self.assertEqual(H.nnz(), 2)                                  # diagonal (objective) only

if __name__ == '__main__':
    unittest.main()   
