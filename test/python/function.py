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

scipy_interpolate = False
try:
  import scipy.interpolate
  scipy.interpolate.RectBivariateSpline
  scipy_interpolate = True
except:
  pass

class Functiontests(casadiTestCase):

  def test_call_empty(self):
    x = SX.sym("x",2)
    fsx = Function("fsx", [x,[]],[x])
    x = MX.sym("x",2)
    fmx1 = Function("fmx1", [x,MX()],[x])
    fmx2 = Function("fmx2", [x,[]],[x])

    for f in [fsx,fmx1,fmx2]:
      f(0,0)

      X = MX.sym("X",2)
      F = f(X,MX())
      g = Function("g", [X],[F])

      g(0)

    x = SX.sym("x",2)
    fsx = Function("fsx", [x],[x,[]])
    x = MX.sym("x",2)
    fmx1 = Function("fmx1", [x],[x,MX()])
    fmx2 = Function("fmx2", [x],[x,[]])

    for f in [fsx,fmx1,]:
      f(0)

      X = MX.sym("X",2)
      F = f(X)
      g = Function("g", [X],F)

      g(0)

  def test_MX_funSeed(self):
    self.message("MX_funSeed")
    x1 = MX.sym("x",2)
    y1 = MX.sym("y")
    x2 = MX.sym("x",2)
    y2 = MX.sym("y")
    p= Function("p", [x1,y1,x2,y2],[sin(x1) + y1,sin(x2) + y2])

    n1 = DM([4,5])
    N1 = 3
    n2 = DM([5,7])
    N2 = 8

    out = p(n1,N1,n2,N2)

    self.checkarray(sin(n1)+N1,out[0],"output")
    self.checkarray(sin(n2)+N2,out[1],"output")
  def test_segfault(self):
    f = Function()
    with self.assertRaises(Exception):
      f.stats()
  def test_issue304(self):
    self.message("regression test for #304") # this code used to segfault
    x = SX.sym("x")

    f = Function("f", [x],[x**2,x**3])

    X = MX.sym("X")

    z=f(X)

    g = Function("g", [X], z).expand()

  def test_jacobian(self):
    x = SX.sym("x",3,1)
    y = SX.sym("y",2,1)

    f = Function("f", [x,y],[x**2,y,x*y[0]])

    g = f.jacobian_old(0, 0)

    self.assertEqual(g.n_in(),f.n_in())
    self.assertEqual(g.n_out(),f.n_out()+1)

  def test_xfunction(self):
    x = SX.sym("x",3,1)
    y = SX.sym("y",2,1)

    f = Function("f", [x,y],[x**2,y,x*y[0]])

    X = MX.sym("x",3,1)
    Y = MX.sym("y",2,1)

    F = Function("F", [X,Y],[X**2,Y,X*Y[0]])

    self.checkfunction(f,F,inputs=[[0.1,0.7,1.3],[7.1,2.9]],sens_der=False,evals=False)


  @memory_heavy()
  def test_jacobians(self):

    x = SX.sym("x")

    self.assertEqual(jacobian(5,x).nnz(),0)


    def test(sp):
      x = SX.sym("x",sp.size2())
      sp2 = jacobian(mtimes(DM.ones(sp),x),x).sparsity()
      self.checkarray(sp.row(),sp2.row());
      self.checkarray(sp.colind(),sp2.colind());

    for i in range(5):
      test(Sparsity.lower(i))
      test(Sparsity.lower(i).T)
      test(Sparsity.dense(i,i))
      test(Sparsity.diag(i))

    for i in [63,64,65,127,128,129]:
      d = Sparsity.diag(i)
      test(d)

      test(d + Sparsity.rowcol([0],[5],i,i))

      b = Sparsity.band(i,-1) + Sparsity.band(i,1)
      test(b + Sparsity.rowcol([0],[5],i,i))

    m = IM.ones(Sparsity.diag(129))
    m[:50,0] = 1
    m[60:,0] = 1
    m[6:9,6] = 1
    m[9,9:12] = 1

    sp = m[:,:120].sparsity()

    test(sp)
    #test(sp.T)

    m = IM.ones(Sparsity.diag(64))
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
      test(Sparsity.lower(i))
      test(Sparsity.lower(i).T)

    for n in ([63,64,65,127,128,129] if args.run_slow else [63,64,65]):
      for m in ([63,64,65,127,128,129] if args.run_slow else [63,64,65]):
        print((n,m))
        sp = Sparsity.dense(n,m)

        test(sp)

        random.seed(0)

        I = IM.ones(sp)
        for i in range(n):
          for j in range(m):
            if random.random()<0.5:
              I[i,j] = 0
        I = sparsify(I)

        sp_holes = I.sparsity()

        test(sp_holes)

        z = IM(sp_holes.size1(), sp_holes.size2())

        R = 5
        v = []
        for r in range(R):
          h = [z]*5
          h[r] = I
          v.append(horzcat(*h))
        d = vertcat(*v)

        test(d.sparsity())

  @memory_heavy()
  def test_hessians(self):
    def test(sp):
      x = SX.sym("x",sp.size2())
      self.assertTrue(sp==sp.T)
      f = Function("f", [x],[mtimes([x.T,DM.ones(sp),x])])
      J = f.hessian_old(0, 0)
      sp2 = J.sparsity_out(0)
      self.checkarray(sp.row(),sp2.row())
      self.checkarray(sp.colind(),sp2.colind())

    A = IM([[1,1,0,0,0,0],[1,1,1,0,1,1],[0,1,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,1,0,1,0,1]])
    A = sparsify(A)
    C = A.sparsity()

    test(C)

    A = IM([[1,0,0,0,0,0],[0,1,1,0,1,1],[0,1,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,1,0,1,0,1]])
    A = sparsify(A)
    C = A.sparsity()

    test(C)

    A = IM([[1,0,0,0,0,0],[0,1,0,0,1,1],[0,0,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,1,0,1,0,1]])
    A = sparsify(A)
    C = A.sparsity()

    test(C)

    A = IM([[0,0,0,0,0,0],[0,1,0,0,1,1],[0,0,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,1,0,1,0,1]])
    A = sparsify(A)
    C = A.sparsity()

    test(C)

    A = IM([[0,0,0,0,0,0],[0,1,0,0,1,0],[0,0,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,0,0,1,0,1]])
    A = sparsify(A)
    C = A.sparsity()

    test(C)


    for i in [63,64,65,100,127,128,129]:
      d = Sparsity.diag(i)
      test(d)

      test(d + Sparsity.rowcol([0],[5],i,i) + Sparsity.rowcol([5],[0],i,i))

      b = Sparsity.band(i,-1) + Sparsity.band(i,1)
      test(b)
      test(b + Sparsity.rowcol([0],[5],i,i) + Sparsity.rowcol([5],[0],i,i))

      d = Sparsity.dense(i,i)
      test(d)

      d = Sparsity.diag(i) + Sparsity.triplet(i,i,[0]*i,list(range(i)))+Sparsity.triplet(i,i,list(range(i)),[0]*i)
      test(d)


      sp = Sparsity.dense(i,i)

      random.seed(0)

      I = IM.ones(sp)
      for ii in range(i):
        for jj in range(i):
          if random.random()<0.5:
            I[ii,jj] = 0
            I[jj,ii] = 0
      I = sparsify(I)

      sp_holes = I.sparsity()

      test(sp_holes)

      z = IM(sp_holes.size1(), sp_holes.size2())

      R = 5
      v = []
      for r in range(R):
        h = [z]*5
        h[r] = I
        v.append(horzcat(*h))
      d = vertcat(*v)

      test(d.sparsity())

  def test_customIO(self):
    x = SX.sym("x")
    f = Function('f',[x],[x*x, x],["i0"], ["foo","bar"])

    ret = f(i0=12)

    self.checkarray(DM([144]),ret["foo"])
    self.checkarray(DM([12]),ret["bar"])


    with self.assertRaises(Exception):
      f_out["baz"]

    ret = f(i0=SX(12))
    self.checkarray(ret["foo"],DM([144]))
    self.checkarray(ret["bar"],DM([12]))
    with self.assertRaises(Exception):
      self.checkarray(ret["baz"],DM([12]))

  def test_derivative_simplifications(self):

    n = 1
    x = SX.sym("x",n)

    M = Function("M", [x],[mtimes((x-DM(list(range(n)))),x.T)])

    P = MX.sym("P",n,n)
    X = MX.sym("X",n)

    M_X= M(X)

    Pf = Function("P", [X, P], [mtimes(M_X,P)])

    P_P = Pf.jacobian_old(1, 0)

    self.assertFalse("derivative" in str(P_P))

  def test_issue1464(self):
    n = 6
    x = SX.sym("x",n)
    u = SX.sym("u")


    N = 9

    rk4 = Function("f",[x,u],[x+u])

    for XX,XFunction in [(SX,Function),(MX,Function)]:

      g = []
      g2 = []


      V = XX.sym("V",(N+1)*n+N)
      VX,VU = vertsplit(V,[0,(N+1)*n,(N+1)*n+N])

      VXk = vertsplit(VX,n)
      VUk = vertsplit(VU,1)

      for k in range(N):

          xf = rk4(VXk[k],VUk[k])

          xfp = vertsplit(xf,int(n/2))
          vp = vertsplit(VXk[k+1],int(n/2))

          g.append(xfp[0] - vp[0])
          g.append(xfp[1] - vp[1])

          g2.append(xf-VXk[k+1])

      for i in range(2):
        f = XFunction("nlp",[V],[vertcat(*g)],{"ad_weight_sp":i})

        assert f.sparsity_jac(0, 0).nnz()==162

        f2 = XFunction("nlp",[V],[vertcat(*g2)],{"ad_weight_sp":i})

        assert f2.sparsity_jac(0, 0).nnz()==162

  def test_callback(self):
    class mycallback(Callback):
      def __init__(self, name, opts={}):
        Callback.__init__(self)
        self.construct(name, opts)

      def eval(self,argin):
        return [argin[0]**2]

    foo = mycallback("my_f")

    x = MX.sym('x')
    y = foo(x)

    f = Function("f",[x],[y])

    out = f(5)

    self.checkarray(out,25)

  def test_callback_errors(self):
    class mycallback(Callback):
      def __init__(self, name, opts={}):
        Callback.__init__(self)
        self.construct(name, opts)
      def eval(self,argin):
        raise Exception("foobar")

    foo = mycallback("my_f")

    x = MX.sym('x')
    y = foo(x)

    f = Function("f",[x],[y])

    try:
      f(3)
    except Exception as e:
      self.assertTrue("foobar" in str(e))

  def test_mapdict(self):
    x = SX.sym("x")
    y = SX.sym("y",2)
    z = SX.sym("z",2,2)
    v = SX.sym("z",Sparsity.upper(3))

    fun = Function("f",{"x":x,"y":y,"z":z,"v":v,"I":mtimes(z,y)+x,"II":sin(y*x).T,"III":v/x},["x","y","z","v"],["I","II","III"])

    n = 2

    X = [MX.sym("x") for i in range(n)]
    Y = [MX.sym("y",2) for i in range(n)]
    Z = [MX.sym("z",2,2) for i in range(n)]
    V = [MX.sym("z",Sparsity.upper(3)) for i in range(n)]

    res = fun.map(n).call({"x":horzcat(*X),"y":horzcat(*Y),"z":horzcat(*Z),"v":horzcat(*V)})

    res2 = fun.map(n).call([horzcat(*X),horzcat(*Y),horzcat(*Z),horzcat(*V)])

    F = Function("F",X+Y+Z+V,res2)
    F2 = Function("F",X+Y+Z+V,[res["I"],res["II"],res["III"]])

    np.random.seed(0)
    X_ = [ DM(i.sparsity(),np.random.random(i.nnz())) for i in X ]
    Y_ = [ DM(i.sparsity(),np.random.random(i.nnz())) for i in Y ]
    Z_ = [ DM(i.sparsity(),np.random.random(i.nnz())) for i in Z ]
    V_ = [ DM(i.sparsity(),np.random.random(i.nnz())) for i in V ]

    self.checkfunction(F,F2,inputs=X_+Y_+Z_+V_,jacobian=False,hessian=False,evals=False)

  @memory_heavy()
  def test_map_node(self):
    x = SX.sym("x")
    y = SX.sym("y",2)
    z = SX.sym("z",2,2)
    v = SX.sym("z",Sparsity.upper(3))

    fun = Function("f",[x,y,z,v],[mtimes(z,y)+x,sin(y*x).T,v/x])

    n = 2

    X = [MX.sym("x") for i in range(n)]
    Y = [MX.sym("y",2) for i in range(n)]
    Z = [MX.sym("z",2,2) for i in range(n)]
    V = [MX.sym("z",Sparsity.upper(3)) for i in range(n)]

    for parallelization in ["serial","openmp","unroll","inline","thread"] if args.run_slow else ["serial"]:
        print(parallelization)
        res = fun.map(n, parallelization).call([horzcat(*x) for x in [X,Y,Z,V]])


        F = Function("F",X+Y+Z+V,list(map(sin,res)))

        resref = [[] for i in range(fun.n_out())]
        for r in zip(X,Y,Z,V):
          for i,e in enumerate(map(sin,fun.call(r))):
            resref[i] = resref[i] + [e]

        Fref = Function("F",X+Y+Z+V,[horzcat(*x) for x in resref])

        np.random.seed(0)
        X_ = [ DM(i.sparsity(),np.random.random(i.nnz())) for i in X ]
        Y_ = [ DM(i.sparsity(),np.random.random(i.nnz())) for i in Y ]
        Z_ = [ DM(i.sparsity(),np.random.random(i.nnz())) for i in Z ]
        V_ = [ DM(i.sparsity(),np.random.random(i.nnz())) for i in V ]

        for f in [F, F.expand('expand_'+F.name())]:

          self.checkfunction(f,Fref,inputs=X_+Y_+Z_+V_,sparsity_mod=args.run_slow)

  def test_map_node_light(self):
    x = SX.sym("x")
    y = SX.sym("y",2)
    z = SX.sym("z",2,2)
    v = SX.sym("z",Sparsity.upper(3))

    fun = Function("f",[x,y,z,v],[mtimes(z,y)+x,sin(y*x).T,v/x])

    n = 2

    X = [MX.sym("x") for i in range(n)]
    Y = [MX.sym("y",2) for i in range(n)]
    Z = [MX.sym("z",2,2) for i in range(n)]
    V = [MX.sym("z",Sparsity.upper(3)) for i in range(n)]

    for parallelization in ["serial","openmp","unroll","inline","thread"]:
        print(parallelization)
        res = fun.map(n, parallelization).call([horzcat(*x) for x in [X,Y,Z,V]])


        F = Function("F",X+Y+Z+V,list(map(sin,res)))

        resref = [[] for i in range(fun.n_out())]
        for r in zip(X,Y,Z,V):
          for i,e in enumerate(map(sin,fun.call(r))):
            resref[i] = resref[i] + [e]

        Fref = Function("F",X+Y+Z+V,[horzcat(*x) for x in resref])

        np.random.seed(0)
        X_ = [ DM(i.sparsity(),np.random.random(i.nnz())) for i in X ]
        Y_ = [ DM(i.sparsity(),np.random.random(i.nnz())) for i in Y ]
        Z_ = [ DM(i.sparsity(),np.random.random(i.nnz())) for i in Z ]
        V_ = [ DM(i.sparsity(),np.random.random(i.nnz())) for i in V ]

        for f in [F, F.expand('expand_'+F.name())]:
          self.checkfunction_light(f,Fref,inputs=X_+Y_+Z_+V_,)

  def test_map_node_n_threads(self):
    x = SX.sym("x")
    y = SX.sym("y",2)
    z = SX.sym("z",2,2)
    v = SX.sym("z",Sparsity.upper(3))

    fun = Function("f",[x,y,z,v],[mtimes(z,y)+x,sin(y*x).T,v/x])

    X_ = [ DM(x.sparsity(),np.random.random(x.nnz())) for i in range(10) ]
    Y_ = [ DM(y.sparsity(),np.random.random(y.nnz())) for i in range(10) ]
    Z_ = [ DM(z.sparsity(),np.random.random(z.nnz())) for i in range(10) ]
    V_ = [ DM(v.sparsity(),np.random.random(v.nnz())) for i in range(10) ]


    print(fun.map(3,"thread",2))


    self.checkfunction_light(fun.map(2,"thread",1),fun.map(2),inputs=[hcat(X_[:2]),hcat(Y_[:2]),hcat(Z_[:2]),hcat(V_[:2])])
    self.checkfunction_light(fun.map(3,"thread",1),fun.map(3),inputs=[hcat(X_[:3]),hcat(Y_[:3]),hcat(Z_[:3]),hcat(V_[:3])])
    self.checkfunction_light(fun.map(3,"thread",2),fun.map(3),inputs=[hcat(X_[:3]),hcat(Y_[:3]),hcat(Z_[:3]),hcat(V_[:3])])
    self.checkfunction_light(fun.map(4,"thread",2),fun.map(4),inputs=[hcat(X_[:4]),hcat(Y_[:4]),hcat(Z_[:4]),hcat(V_[:4])])
    self.checkfunction_light(fun.map(4,"thread",5),fun.map(4),inputs=[hcat(X_[:4]),hcat(Y_[:4]),hcat(Z_[:4]),hcat(V_[:4])])

  @memory_heavy()
  def test_mapsum(self):
    x = SX.sym("x")
    y = SX.sym("y",2)
    z = SX.sym("z",2,2)
    v = SX.sym("z",Sparsity.upper(3))

    fun = Function("f",[x,y,z,v],[mtimes(z,y)+x,sin(y*x).T,v/x])

    n = 2

    X = [MX.sym("x") for i in range(n)]
    Y = [MX.sym("y",2) for i in range(n)]
    Z = [MX.sym("z",2,2) for i in range(n)]
    V = [MX.sym("z",Sparsity.upper(3)) for i in range(n)]

    zi = 0
    for Z_alt in [Z,[MX()]*3]:
      zi+= 1
      for parallelization in ["serial","openmp","unroll","thread"]:
        res = fun.mapsum([horzcat(*x) for x in [X,Y,Z_alt,V]],parallelization) # Joris - clean alternative for this?

        for ad_weight_sp in [0,1]:
          F = Function("F",X+Y+Z+V,list(map(sin,res)),{"ad_weight": 0,"ad_weight_sp":ad_weight_sp})

          resref = [0 for i in range(fun.n_out())]
          for r in zip(X,Y,Z_alt,V):
            for i,e in enumerate(fun.call(r)):
              resref[i] = resref[i] + e

          Fref = Function("F",X+Y+Z+V,list(map(sin,resref)))

          np.random.seed(0)
          X_ = [ DM(i.sparsity(),np.random.random(i.nnz())) for i in X ]
          Y_ = [ DM(i.sparsity(),np.random.random(i.nnz())) for i in Y ]
          Z_ = [ DM(i.sparsity(),np.random.random(i.nnz())) for i in Z ]
          V_ = [ DM(i.sparsity(),np.random.random(i.nnz())) for i in V ]

          inputs = X_+Y_+Z_+V_

          if parallelization!="thread":
            self.check_codegen(F,inputs=inputs)

          for f in [F,toSX_fun(F)]:
            self.checkfunction(f,Fref,inputs=inputs,sparsity_mod=args.run_slow)


  @memory_heavy()
  def test_mapsum2(self):
    x = SX.sym("x")
    y = SX.sym("y",2)
    z = SX.sym("z",2,2)
    v = SX.sym("z",Sparsity.upper(3))

    fun = Function("f",[x,y,z,v],[mtimes(z,y)+x,sin(y*x).T,v/x])

    n = 2

    X = [MX.sym("x") for i in range(n)]
    Y = [MX.sym("y",2) for i in range(n)]
    Z = MX.sym("z",2,2)
    V = MX.sym("z",Sparsity.upper(3))

    for Z_alt in [Z]:

      for parallelization in ["serial","openmp","unroll","thread"]:

        for ad_weight_sp in [0,1]:
          for ad_weight in [0,1]:
            F = fun.map("map",parallelization,n,[2,3],[0],{"ad_weight_sp":ad_weight_sp,"ad_weight":ad_weight})

            resref = [0 for i in range(fun.n_out())]
            acc = 0
            bl = []
            cl = []
            for r in zip(X,Y,[Z_alt]*n,[V]*n):
              a,b,c= fun(*r)
              acc = acc + a
              bl.append(b)
              cl.append(c)

            Fref = Function("F",[horzcat(*X),horzcat(*Y),Z,V],[acc,horzcat(*bl),horzcat(*cl)])

            np.random.seed(0)
            X_ = [ DM(i.sparsity(),np.random.random(i.nnz())) for i in X ]
            Y_ = [ DM(i.sparsity(),np.random.random(i.nnz())) for i in Y ]
            Z_ = DM(Z.sparsity(),np.random.random(Z.nnz()))
            V_ = DM(V.sparsity(),np.random.random(V.nnz()))

            inputs = [horzcat(*X_),horzcat(*Y_),Z_,V_]

            if parallelization!="thread":
              self.check_codegen(F,inputs=inputs)

            for f in [F,toSX_fun(F)]:
              self.checkfunction(f,Fref,inputs=inputs,sparsity_mod=args.run_slow)

  def test_repmatnode(self):
    x = MX.sym("x",2)

    y = sin(repmat(x**2,1,3))

    z = MX.sym("y",2,2)

    F = Function("f",[x,z],[sum2(sum1(y))])

    x = SX.sym("x",2)

    y = sin(repmat(x**2,1,3))
    z = SX.sym("y",2,2)

    Fref = Function("f",[x,z],[sum2(sum1(y))])

    x0 = DM([1,7])
    x1 = DM([[3,0],[2,4]])

    self.check_codegen(F,inputs=[x0,x1])
    self.checkfunction(F,Fref,inputs=[x0,x1])

  def test_repsumnode(self):

    x = MX.sym("x",2)
    z = MX.sym("y",2,2)

    F = Function("f",[x,z],[sin(repsum((x**2).T,1,2)),(cos(x**2)*2*x).T])

    x = SX.sym("x",2)
    z = SX.sym("y",2,2)


    Fref = Function("f",[x,z],[sin(repsum((x**2).T,1,2)),(cos(x**2)*2*x).T])

    x0 = DM([1,7])
    x1 = DM([[3,0],[2,4]])

    self.check_codegen(F,inputs=[x0,x1])
    self.checkfunction(F,Fref,inputs=[x0,x1])

  def test_unknown_options(self):
    x = SX.sym("x")

    with self.assertRaises(Exception):
      f = SXFunction("f", [x],[x],{"fooo": False})

    with self.assertRaises(Exception):
      f = SXFunction("f", [x],[x],{"ad_weight": "foo"})

    if not has_nlpsol("ipopt"):
      return

  @known_bug()
  def test_unknown_options_stringvector(self):
    x = SX.sym("x")
    solver = nlpsol("mysolver", "ipopt", {"x":x,"f":x**2}, {"monitor": ["eval_f"]})
    with self.assertRaises(Exception):
      solver = nlpsol("mysolver", "ipopt", {"x":x,"f":x**2}, {"monitor": ["abc"]})

  @memory_heavy()
  def test_mapaccum(self):

    x = SX.sym("x",2)
    y = SX.sym("y")
    z = SX.sym("z",2,2)
    v = SX.sym("v",Sparsity.upper(3))

    fun = Function("f",[x,y,z,v],[mtimes(z,x)+y,sin(y*x).T,v/y])

    n = 2

    X = MX.sym("x",x.sparsity())
    Y = [MX.sym("y",y.sparsity()) for i in range(n)]
    Z = [MX.sym("z",z.sparsity()) for i in range(n)]
    V = [MX.sym("v",v.sparsity()) for i in range(n)]

    np.random.seed(0)
    X_ = DM(x.sparsity(),np.random.random(x.nnz()))
    Y_ = [ DM(i.sparsity(),np.random.random(i.nnz())) for i in Y ]
    Z_ = [ DM(i.sparsity(),np.random.random(i.nnz())) for i in Z ]
    V_ = [ DM(i.sparsity(),np.random.random(i.nnz())) for i in V ]

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
        Fref = Function("f",[X,horzcat(*Y),horzcat(*Z),horzcat(*V)],[horzcat(*Xps),horzcat(*Y0s),horzcat(*Y1s)])
        inputs = [X_,horzcat(*Y_),horzcat(*Z_),horzcat(*V_)]

        for f in [F,toSX_fun(F)]:

          self.checkfunction(f,Fref,inputs=inputs)
          self.check_codegen(f,inputs=inputs)

    fun = Function("f",[y,x,z,v],[mtimes(z,x)+y+c.trace(v)**2,sin(y*x).T,v/y])

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

        Fref = Function("f",[horzcat(*Y),X,horzcat(*Z),V[0]],[horzcat(*Xps),horzcat(*Y0s),horzcat(*Vps)])
        inputs = [horzcat(*Y_),X_,horzcat(*Z_),V_[0]]

        for f in [F,toSX_fun(F)]:
          self.checkfunction(f,Fref,inputs=inputs)
          self.check_codegen(f,inputs=inputs)

  def test_mapaccum_schemes(self):

    x = SX.sym("x",2)
    y = SX.sym("y")
    z = SX.sym("z",2,2)
    v = SX.sym("v",Sparsity.upper(3))

    fun = Function("f",[y,z,x,v],[mtimes(z,x)+y,sin(y*x).T,v/y],["y","z","x","v"],["out0","out1","out2"])

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

    fun = Function("f",[x,y,z,v],[mtimes(z,x)+y,sin(y*x).T,v/y],["x","y","z","v"],["out0","out1","out2"])

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
    x = SX.sym("x")
    y = x**2
    try:
        depends_on(x,y)
    except Exception as e:
        s = str(e)
    self.assertTrue("not symbolic" in s)
    try:
        Function("f",[y],[x])
    except Exception as e:
        s = str(e)
    self.assertTrue("not symbolic" in s)

  def test_1d_interpolant(self):
    grid = [[0, 1, 1.5, 2, 3]]
    values = [0, 1, 2, 5, 3]
    F = interpolant('F', 'linear', grid, values)
    def same(a, b): return abs(float(a)-b)<1e-8
    pairs = [
      (3.4,3-0.4*2),
      (2.4,5-0.4*2),
      (1.6,2+3*0.1/0.5),
      (1.4,1+0.4/0.5),
      (0.4,0.4),
      (-.6,-0.6)
    ]

    X = MX.sym("x")

    J = Function("F",[X],[F(X)])

    for a,r in pairs:
      self.assertTrue(same(F(a), r))
      self.check_codegen(F,inputs=[a])


    X = MX.sym("x")

    J = Function("F",[X],[jacobian(F(X),X)])

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

  def test_2d_interpolant(self):
    grid = [[0, 1, 4, 5],
            [0, 2, 3]]

    values = [0,   1,  8,  3,
              10, -11, 12, 13,
              20, 31, -42, 53]
    F = interpolant('F', 'linear', grid, values)


    a0 = -11+0.4*(31+11)
    a1 = 12+0.4*(-42-12)
    pairs = [
      (vertcat(1,2), -11),
      (vertcat(1,3), 31),
      (vertcat(4,2), 12),
      (vertcat(4,3), -42),

      (vertcat(1,2.4), a0),
      (vertcat(4,2.4), a1),

      (vertcat(3,2), -11+2.0/3*(12+11)),
      (vertcat(3,3), 31+2.0/3*(-42-31)),

      (vertcat(3,2.4), a0+2.0/3*(a1-a0))
    ]

    for a,r in pairs:
      self.checkarray(F(a), r)
      self.check_codegen(F,inputs=[a])


    X = MX.sym("x",2)

    J = Function("J",[X],[jacobian(F(X),X)])

    jx0 = (12+11)/3.0
    jx1 = (-42-31)/3.0
    jx2 = (13-12)
    jx3 = (53+42)

    jy0 = 31+11
    jy1 = -42-12

    pairs = [
      (vertcat(1,2), vertcat(jx0,jy0)),
      (vertcat(1,3), vertcat(jx1,jy0)),
      (vertcat(4,2), vertcat(jx2,jy1)),
      (vertcat(4,3), vertcat(jx3,jy1)),

      (vertcat(1,2.4), vertcat(jx0+(jx1-jx0)*0.4, 31+11)),
      (vertcat(4,2.4), vertcat(jx2+(jx3-jx2)*0.4, -42-12)),

      (vertcat(3,2), vertcat(jx0,jy0+(jy1-jy0)*2.0/3)),
      (vertcat(3,3), vertcat(jx1,jy0+(jy1-jy0)*2.0/3)),

      (vertcat(3,2.4), vertcat(jx0+(jx1-jx0)*0.4,jy0+(jy1-jy0)*2.0/3)),

    ]

    for a,r in pairs:
      self.checkarray(J(a).T, r)
      self.check_codegen(J,inputs=[a])

  def test_1d_interpolant_uniform(self):
    grid = [[0, 1, 2]]
    values = [0, 1, 2]
    for opts in [{"lookup_mode": ["linear"]},{"lookup_mode": ["exact"]},{"lookup_mode": ["binary"]}]:
      F = interpolant('F', 'linear', grid, values, opts)
      def same(a, b): return abs(float(a)-b)<1e-8
      self.assertTrue(same(F(2.4), 2.4))
      self.assertTrue(same(F(1.4), 1.4))
      self.assertTrue(same(F(0), 0))
      self.assertTrue(same(F(1), 1))
      self.assertTrue(same(F(2), 2))
      self.assertTrue(same(F(6), 6))
      self.assertTrue(same(F(0.4), 0.4))
      self.assertTrue(same(F(-.6), -.6))

      F = interpolant('F', 'linear', [np.linspace(0,1,7)], list(range(7)), {"lookup_mode": ["exact"]})

    grid = [[2, 4, 6]]
    values = [10, 7, 1]
    for opts in [{"lookup_mode": ["linear"]},{"lookup_mode": ["exact"]},{"lookup_mode": ["binary"]}]:
      F = interpolant('F', 'linear', grid, values, opts)
      def same(a, b): return abs(float(a)-b)<1e-8
      self.assertTrue(same(F(1), 11.5))
      self.assertTrue(same(F(2), 10))
      self.assertTrue(same(F(3), 8.5))
      self.assertTrue(same(F(4), 7))
      self.assertTrue(same(F(5), 4))
      self.assertTrue(same(F(6), 1))
      self.assertTrue(same(F(7), -2))

      F = interpolant('F', 'linear', [np.linspace(0,1,7)], list(range(7)), {"lookup_mode": ["exact"]})

  def test_2d_interpolant_uniform(self):
    grid = [[0, 1, 2], [0, 1, 2]]
    values = [0, 1, 2, 10, 11, 12, 20, 21, 22]
    for opts in [{"lookup_mode": ["linear","linear"]},{"lookup_mode": ["exact","exact"]},{"lookup_mode": ["binary","binary"]}]:
      F = interpolant('F', 'linear', grid, values, opts)
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
  def test_2d_bspline(self):
    import scipy.interpolate
    np.random.seed(0)

    d_knots = [list(np.linspace(0,1,5)),list(np.linspace(0,1,6))]

    data = np.random.random([len(e) for e in d_knots])
    r = np.meshgrid(*d_knots,indexing='ij')

    xyz = np.vstack(e.ravel(order='F') for e in r).ravel(order='F')

    d_flat = data.ravel(order='F')

    LUT = casadi.interpolant('name','bspline',d_knots,d_flat)
    LUTJ = LUT.jacobian_old(0, 0)
    LUTH = LUT.hessian_old(0, 0)

    self.check_codegen(LUT, [vertcat(0.2,0.3)])
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
          r = blockcat([[interp.ev(x,y, 2, 0),interp.ev(x,y, 1, 1)],[interp.ev(x,y, 1, 1), interp.ev(x,y, 0, 2)]])
        except:
          r = None
        if r is not None:
          self.checkarray(m,r)

  @skip(not scipy_interpolate)
  def test_1d_bspline(self):
    import scipy.interpolate
    np.random.seed(0)

    d_knots = [list(np.linspace(0,1,5))]

    data = np.random.random([len(e) for e in d_knots])
    r = np.array(d_knots)

    xyz = np.vstack(e.ravel(order='F') for e in r).ravel(order='F')

    d_flat = data.ravel(order='F')

    LUT = casadi.interpolant('name','bspline',d_knots,d_flat)
    self.check_codegen(LUT, [0.2])
    LUTJ = LUT.jacobian_old(0, 0)
    LUTH = LUT.hessian_old(0, 0)

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
    x = MX.sym("x")
    y = MX.sym("y")

    num_inputs = [0.2,0.7]

    g = Function("g", [x,y],[sin(x+3*y)])

    class Fun(Callback):
        # sin(x+3*y)

        def __init__(self):
          Callback.__init__(self)
          self.construct("Fun", {})
        def get_n_in(self): return 2
        def get_n_out(self): return 1

        def eval(self,arg):
          x = arg[0]
          y = arg[1]
          z0 = 3*y
          z1 = x+z0
          z2 = sin(z1)
          return [z2]

        def has_forward(self,nfwd): return False
        def has_reverse(self,nadj): return False

        def has_jacobian(self): return True

        def get_jacobian(self, name, inames, onames, opts):
          x = SX.sym("x")
          y = SX.sym("y")
          out_g = SX.sym('out_g', Sparsity(1,1))
          J = Function(name, [x,y,out_g],[horzcat(cos(x+3*y),3*cos(x+3*y))], inames, onames, opts)
          return J

    f = Fun()

    self.checkfunction(f,g,inputs=num_inputs,fwd=False,adj=False,indirect=False)


  def test_Callback_errors(self):

    class Fun(Callback):

        def __init__(self):
          Callback.__init__(self)
          self.construct("Fun", {})
        def get_n_in(self): return 2
        def get_n_out(self): return 1

        def get_sparsity_in(i):
          return 4

        def eval(self,arg):
          x = arg[0]
          y = arg[1]

          z0 = 3*y
          z1 = x+z0
          z2 = sin(z1)
          return [z2]

    try:
      f = Fun()
    except Exception as e:
      s = str(e)
      print(s)
    self.assertTrue("get_sparsity_in" in s)

  def test_Callback(self):

    x = MX.sym("x")
    y = MX.sym("y")

    num_inputs = [0.2,0.7]

    g = Function("g", [x,y],[sin(x+3*y)])

    # Simple syntax
    def getP(indirect=True):

      class Fun(Callback):

        def __init__(self):
          Callback.__init__(self)
          self.construct("Fun", {})

        def get_n_in(self): return 2
        def get_n_out(self): return 1

        def eval(self,arg):
          x = arg[0]
          y = arg[1]

          z0 = 3*y
          z1 = x+z0
          z2 = sin(z1)
          return [z2]

      self.cb = Fun()

      if not indirect:
        return self.cb

      f = Function("f", [x,y],[self.cb(x,y)])

      return f

    for indirect in [True,False]:
      f = getP(indirect=indirect)
      self.checkfunction(f,g,inputs=num_inputs,sens_der=False,jacobian=False,gradient=False,hessian=False,evals=False)

      with self.assertRaises(Exception):
        f.gradient()

      with self.assertRaises(Exception):
        f.jacobian_old(0, 0)

      with self.assertRaises(Exception):
        f.forward(1)
      with self.assertRaises(Exception):
        f.reverse(1)

  def test_Callback_dimcheck(self):
      class Fun(Callback):
        def __init__(self):
          Callback.__init__(self)
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
      class Fun(Callback):
        def __init__(self):
          Callback.__init__(self)
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
      class Fun(Callback):
        def __init__(self):
          Callback.__init__(self)
          self.construct("Fun")
        def get_n_in(self): return 2
        def get_n_out(self): return 1

        def eval(self,arg):
          return [DM.zeros(2,2)]
      f = Fun()
      try:
        f(2,3)
      except Exception as e:
        s = str(e)
      self.assertTrue("Shape mismatch" in s)

  def test_Callback_sens(self):
    x = MX.sym("x")
    y = MX.sym("y")

    num_inputs = [0.2,0.7]

    g = Function("g", [x,y],[sin(x+3*y)])

    def getP(has_fwd=True,has_adj=True,indirect=True):

      class Fun(Callback):
        # sin(x+3*y)

        def __init__(self,opts):
          Callback.__init__(self)
          self.construct("Fun", opts)
        def get_n_in(self): return 2
        def get_n_out(self): return 1

        def eval(self,arg):
          x = arg[0]
          y = arg[1]
          z0 = 3*y
          z1 = x+z0
          z2 = sin(z1)
          return [z2]

        if has_fwd:
          def has_forward(self,nfwd): return nfwd==1
          def get_forward(self,nfwd,name,inames,onames,opts):
            assert(nfwd==1)
            class ForwardFun(Callback):
              # sin(x+3*y)

              def __init__(self):
                Callback.__init__(self)
                self.construct(name, {"verbose":True})
              def get_n_in(self): return 2+1+2
              def get_n_out(self): return 1

              def eval(self,arg):
                x,y = arg[0],arg[1]
                z = arg[2]
                seeds = arg[3:]

                z0 = 3*y
                z1 = x+z0
                z2 = sin(z1)

                ret = []

                for i in range(3,len(arg),2):
                  dx = arg[i]
                  dy = arg[i+1]
                  dz0 = 3*dy
                  dz1 = dx+dz0
                  dz2 = cos(z1)*dz1
                  ret.append(dz2)

                return ret

            self.cb_fwd = ForwardFun()
            return self.cb_fwd

        if has_adj:
          def has_reverse(self,nadj): return nadj==1
          def get_reverse(self,nadj,name,inames,onames,opts):
            assert(nadj==1)
            class BackwardFun(Callback):
              # sin(x+3*y)

              def __init__(self):
                Callback.__init__(self)
                self.construct(name, {"verbose":True})
              def get_n_in(self): return 2+1+1
              def get_n_out(self): return 2

              def eval(self,arg):
                x,y = arg[0],arg[1]
                z = arg[2]
                seeds = arg[3:]

                z0 = 3*y
                z1 = x+z0
                z2 = sin(z1)

                ret = []

                for i in range(3,len(arg)):
                  z_bar = arg[i]
                  bx = 0
                  by = 0
                  bz1 = 0
                  bz0 = 0

                  bz2 = z_bar
                  bz1 += bz2*cos(z1)
                  bx+= bz1;bz0+= bz1
                  by+= 3*bz0
                  ret.append(bx)
                  ret.append(by)
                return ret

            self.cb_rev = BackwardFun()
            return self.cb_rev

      opts = {"verbose":True}
      self.cb = Fun(opts)
      f = self.cb

      if not indirect:
        return f

      f = Function("f", [x,y],[f(x,y)])

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

  @requires_nlpsol("ipopt")
  def test_common_specific_options(self):

      x = SX.sym("x")

      nlp = {"x": x, "f": x**2}

      with capture_stdout() as out:
        solver = nlpsol("solver","ipopt",nlp)
      self.assertTrue("nlp_f" not in out[0])
      with capture_stdout() as out:
        solver = nlpsol("solver","ipopt",nlp,{"common_options":{"verbose":True}})
      self.assertTrue("nlp_f" in out[0])
      with capture_stdout() as out:
        solver = nlpsol("solver","ipopt",nlp,{"specific_options":{ "nlp_f" : {"verbose":True}}})
      self.assertTrue("nlp_f" in out[0])
      with capture_stdout() as out:
        solver = nlpsol("solver","ipopt",nlp,{"common_options":{"verbose":True},"specific_options":{ "nlp_f" : {"verbose":False}}})
      self.assertTrue("nlp_f" not in out[0])
      with capture_stdout() as out:
        solver = nlpsol("solver","ipopt",nlp,{"common_options":{"verbose":False},"specific_options":{ "nlp_f" : {"verbose":True}}})
      self.assertTrue("nlp_f" in out[0])

      with capture_stdout() as out:
        solver = nlpsol("solver","ipopt",nlp)
      self.assertTrue(len(out[1])==0)
      with capture_stdout() as out:
        solver = nlpsol("solver","ipopt",nlp,{"specific_options":{ "nlp_foo" : {"verbose":True}}})
      self.assertTrue("Ignoring" + out[1])
      self.assertTrue("nlp_g" in out[1])
      with self.assertRaises(Exception):
        solver = nlpsol("solver","ipopt",nlp,{"specific_options":{ "nlp_foo" : 3}})

  @requires_expm("slicot")
  @memory_heavy()
  def test_expm(self):
      eps = 1e-6
      t = MX.sym('t')
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
        x = MX.sym('x',n)
        As = MX.sym('A',n,n)

        dae = {'x':x,'p':vec(As),'ode':mtimes(As,x)}
        intg = integrator('intg','cvodes',dae,{'reltol':1e-14,'abstol':1e-14})

        Intg = intg.map('identity','serial',n,[1],[])

        out = Intg(x0=DM.eye(n),p=vec(As))
        expmF = Function('expm',[As],[out["xf"]])
        return expmF(A)

      A = MX.sym("A",n,n)
      t = MX.sym("t")
      fr = Function('fr',[A,t],[expm(A*t)])
      f = Function('f',[A,t],[casadi.expm(A*t)])

      self.checkfunction(fr,f,inputs=[Anum, 1.1],digits=8)

      fr = Function('fr',[t],[expm(Anum*t)])
      f = Function('f',[t],[casadi.expm_const(Anum,t)])

      self.checkfunction(fr,f,inputs=[1.1],digits=8)

      JA = jacobian(casadi.expm(A*t),A)
      Jt = jacobian(casadi.expm(A*t),t)

      self.assertTrue(JA.nnz()==n**4)
      self.assertTrue(Jt.nnz()==n**2)

      JA = jacobian(casadi.expm_const(A,t),A)
      Jt = jacobian(casadi.expm_const(A,t),t)

      self.assertTrue(JA.nnz()==0)
      self.assertTrue(Jt.nnz()==n**2)

  def test_conditional(self):

    np.random.seed(5)

    x = MX.sym('x',2,2)
    y = MX.sym('y',2,2)

    sp1 = MX.sym('y',Sparsity.lower(2))
    sp2 = MX.sym('z',Sparsity.diag(2))

    f1 = Function("f",[sp2,x],[x**2,x*sp2])
    f2 = Function("f",[sp1,x],[2*x**2,sin(sp1)])
    f3 = Function("f",[sp1,sp2],[sp1*sp2,sp1+sp2])

    F = Function.conditional("test",[f1,f2], f3)
    Fsx = F.expand()

    A = np.random.random((2,2))
    B = np.random.random((2,2))

    for i in range(-1,3):
      self.checkfunction(F,Fsx,inputs = [i,A,B])
      self.check_codegen(F,inputs=[i,A,B])

  def test_max_num_dir(self):
    x = MX.sym("x",10)

    f = Function("ffff",[x],[mtimes(DM.ones(10,10),x)],{"max_num_dir":4,"verbose":True})
    f = f.expand()


    y = f.call([sin(x)],False,True)[0]



    acc = Function('acc',[x],[y])
    acc = acc.mapaccum('accd',5)


    cons = vec(acc(x))

    os.system("callgrind_control -i on")

    e = jacobian(cons,x)
    os.system("callgrind_control -i off")

    g = Function('g',[x],[e])


    c = CodeGenerator('me')
    c.add(g)
    code= c.dump()

    self.assertTrue("fwd4_ffff" in code)

  def test_max_num_dir(self):
    x = MX.sym("x")

    f = Function("ffff",[x],[sin(x)])
    fm = f.mapaccum("mapaccum",100,1,{"base":4})

    c = CodeGenerator('me')
    c.add(fm)
    code= c.dump()

    self.assertTrue("ffff_acc4_acc4_acc4" in code)

  def test_2d_linear_multiout(self):
    np.random.seed(0)

    d_knots = [list(np.linspace(0,1,5)),list(np.linspace(0,1,6))]

    data0 = np.random.random([len(e) for e in d_knots])
    data1 = np.random.random([len(e) for e in d_knots])
    r = np.meshgrid(*d_knots,indexing='ij')

    xyz = np.vstack(e.ravel(order='F') for e in r).ravel(order='F')

    d_flat0 = data0.ravel(order='F')

    LUT0 = casadi.interpolant('name','linear',d_knots,d_flat0)

    d_flat1 = data1.ravel(order='F')

    LUT1 = casadi.interpolant('name','linear',d_knots,d_flat1)

    data = np.vstack((data0.ravel(order='F'),data1.ravel(order='F'))).ravel(order='F')

    d_flat = data.ravel(order='F')


    LUT = casadi.interpolant('name','linear',d_knots,d_flat)


    x = MX.sym("x")
    y = MX.sym("y")


    LUT_sep = Function('f',[x,y],[vertcat(LUT0(vertcat(x,y)),LUT1(vertcat(x,y)))])
    LUT = Function('f',[x,y],[LUT(vertcat(x,y))])

    self.checkfunction(LUT,LUT_sep, inputs=[0.2,0.333])
    self.check_codegen(LUT,inputs=[0.2,0.333])

  def test_2d_bspline_multiout(self):
    np.random.seed(0)

    d_knots = [list(np.linspace(0,1,5)),list(np.linspace(0,1,6))]

    data0 = np.random.random([len(e) for e in d_knots])
    data1 = np.random.random([len(e) for e in d_knots])
    r = np.meshgrid(*d_knots,indexing='ij')

    xyz = np.vstack(e.ravel(order='F') for e in r).ravel(order='F')

    d_flat0 = data0.ravel(order='F')

    LUT0 = casadi.interpolant('name','bspline',d_knots,d_flat0)

    d_flat1 = data1.ravel(order='F')

    LUT1 = casadi.interpolant('name','bspline',d_knots,d_flat1)

    data = np.vstack((data0.ravel(order='F'),data1.ravel(order='F'))).ravel(order='F')

    d_flat = data.ravel(order='F')


    LUT = casadi.interpolant('name','bspline',d_knots,d_flat)


    x = MX.sym("x")
    y = MX.sym("y")


    LUT_sep = Function('f',[x,y],[vertcat(LUT0(vertcat(x,y)),LUT1(vertcat(x,y)))])
    LUT = Function('f',[x,y],[LUT(vertcat(x,y))])

    self.checkfunction(LUT,LUT_sep, inputs=[0.2,0.333])
    self.check_codegen(LUT,inputs=[0.2,0.333])

  def test_codegen_avoid_stack(self):
    x = SX.sym("x",3,3)
    f = Function('f',[x],[det(x)])
    np.random.seed(0)
    self.check_codegen(f,inputs=[np.random.random((3,3))])
    self.check_codegen(f,inputs=[np.random.random((3,3))], opts={"avoid_stack": True})


  def test_sx_serialize(self):
    x = SX.sym("x")
    y = x+3
    z = sin(y)

    f = Function('f',[x],[z])
    fs = Function.deserialize(f.serialize())

    self.checkfunction(f,fs,inputs=[2])

    x = SX.sym("x")
    y = x+3
    z = sin(y)

    f = Function('f',[x],[z,np.nan,-np.inf,np.inf])
    fs = Function.deserialize(f.serialize())
    self.checkfunction(f,fs,inputs=[2])

    x = SX.sym("x")
    y = SX.sym("y", Sparsity.lower(3))
    z = x+y
    z1 = sparsify(vertcat(z[0],0,z[1]))
    z2 = z.T

    f = Function('f',[x,y],[z1,z2,x**2],["x","y"],["a","b","c"])
    fs = Function.deserialize(f.serialize())

    self.assertEqual(fs.name_in(0), "x")
    self.assertEqual(fs.name_out(0), "a")
    self.assertEqual(fs.name(), "f")

    self.checkfunction(f,fs,inputs=[3.7,np.array([[1,0,0],[2,3,0],[4,5,6]])],hessian=False)


    fs = pickle.loads(pickle.dumps(f))
    self.checkfunction(f,fs,inputs=[3.7,np.array([[1,0,0],[2,3,0],[4,5,6]])],hessian=False)

    x = SX.sym("x")
    p = SX.sym("p")

    f = Function('f',[x],[p])

    with self.assertInException("Cannot serialize SXFunction with free parameters"):
      pickle.loads(pickle.dumps(f))


    x = MX.sym("x")
    f = Function('f',[x],[x**2])

    with self.assertInException("'serialize' not defined for MXFunction"):
      pickle.loads(pickle.dumps(f))

  def test_string(self):
    x=MX.sym("x")
    y=MX.sym("y")
    f = Function('f',[x],[],["x"],[])
    self.assertTrue("(x)->()" in str(f))

    f = Function('f',[],[x],[],["y"])
    self.assertTrue("()->(y)" in str(f))

    f = Function('f',[x],[x**2],["x"],["y"])
    self.assertTrue("(x)->(y)" in str(f))

    f = Function('f',[x,y],[x**2,x*y],["x","y"],["w","z"])
    self.assertTrue("(x,y)->(w,z)" in str(f))

  def test_fold(self):
    x = SX.sym("x",2,2)

    f = Function("f",[x],[sin(x)])

    F = f.fold(10)

    x0 = x
    for i in range(10):
      x = f(x)
    Fref = Function("f",[x0],[x])

    self.checkfunction(F,Fref,inputs=[DM([[1,2],[3,7]])])


  @memory_heavy()
  def test_thread_safety(self):
    x = MX.sym('x')
    y = MX.sym('y')
    f = Function('f', [x, y], [x ** 2 - y])
    for rf in ["newton","fast_newton"]:
      finv = rootfinder('finv', rf, f)

      finv_par = finv.map(50,"unroll").map(4, 'thread',4)
      res = finv_par(numpy.ones(200), numpy.linspace(0, 10, 200))
      self.checkarray(norm_inf(res.T-sqrt(numpy.linspace(0, 10, 200))),0, digits=5)

      finv_par = finv.map(200, 'thread',4)
      res = finv_par(numpy.ones(200), numpy.linspace(0, 10, 200))
      self.checkarray(norm_inf(res.T-sqrt(numpy.linspace(0, 10, 200))),0, digits=5)

  def test_mapped_eval(self):
      x = SX.sym('x')
      y = SX.sym('y', 2)
      f = Function('f', [x,y], [sin(x)*y], ['x','y'], ['r'])
      r1 = f(1, DM([3,4]))
      r2 = f(2, DM([3,4]))
      r3 = f(3, DM([3,4]))
      r_all = f(DM([[1,2,3]]), DM([3,4]))
      self.checkarray(r_all, horzcat(r1,r2,r3), "Mapped evaluation (DM)")

      z = MX.sym("z", 1, 3);
      rz = f(z, DM([3,4]))
      F = Function('F', [z], [rz]);
      r_mx = F(DM([[1,2,3]]))
      self.checkarray(r_all, r_mx, "Mapped evaluation (MX)")


if __name__ == '__main__':
    unittest.main()
