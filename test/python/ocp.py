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

class OCPtests(casadiTestCase):
  @requires_nlpsol("ipopt")
  def testdiscrete(self):
    self.message("Linear-quadratic problem, discrete, using IPOPT")
    # inspired by www.cs.umsl.edu/~janikow/publications/1992/GAforOpt/text.pdf
    a=1.0
    b=1.0
    q=1.0
    s=1.0
    r=1.0
    x0=100

    N=100

    X=SX.sym("X",N+1)
    U=SX.sym("U",N)

    V = vertcat(*[X,U])

    cost = 0
    for i in range(N):
      cost = cost + s*X[i]**2+r*U[i]**2
    cost = cost + q*X[N]**2

    nlp = {'x':V, 'f':cost, 'g':vertcat(*[X[0]-x0,X[1:,0]-(a*X[:N,0]+b*U)])}
    opts = {}
    opts["ipopt.tol"] = 1e-5
    opts["ipopt.hessian_approximation"] = "limited-memory"
    opts["ipopt.max_iter"] = 100
    opts["ipopt.print_level"] = 0
    solver = nlpsol("solver", "ipopt", nlp, opts)
    solver_in = {}
    solver_in["lbx"]=[-1000 for i in range(V.nnz())]
    solver_in["ubx"]=[1000 for i in range(V.nnz())]
    solver_in["lbg"]=[0 for i in range(N+1)]
    solver_in["ubg"]=[0 for i in range(N+1)]
    solver_out = solver(**solver_in)
    ocp_sol=solver_out["f"][0]
    # solve the ricatti equation exactly
    K = q+0.0
    for i in range(N):
      K = s+r*a**2*K/(r+b**2*K)
    exact_sol=K * x0**2
    self.assertAlmostEqual(ocp_sol,exact_sol,10,"Linear-quadratic problem solution using IPOPT")

  @requires_nlpsol("ipopt")
  def test_singleshooting(self):
    self.message("Single shooting")
    p0 = 0.2
    y0= 1
    yc0=dy0=0
    te=0.4

    t=SX.sym("t")
    q=SX.sym("y",2,1)
    p=SX.sym("p",1,1)
    # y
    # y'
    dae={'x':q, 'p':p, 't':t, 'ode':vertcat(*[q[1],p[0]+q[1]**2 ])}
    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    opts["verbose"] = False
    opts["steps_per_checkpoint"] = 10000
    opts["t0"] = 0
    opts["tf"] = te
    integrator = casadi.integrator("integrator", "cvodes", dae, opts)

    var = MX.sym("var",2,1)
    par = MX.sym("par",1,1)
    parMX= par

    q0   = vertcat(*[var[0],par])
    par  = var[1]
    qend = integrator(x0=q0, p=par)["xf"]

    parc = MX(0)

    f = Function('f', [var,parMX],[qend[0]])
    nlp = {'x':var, 'f':-f(var,parc)}
    opts = {}
    opts["ipopt.tol"] = 1e-12
    opts["ipopt.hessian_approximation"] = "limited-memory"
    opts["ipopt.max_iter"] = 10
    opts["ipopt.derivative_test"] = "first-order"
    opts["ipopt.print_level"] = 0
    solver = nlpsol("solver", "ipopt", nlp, opts)
    solver_in = {}
    solver_in["lbx"]=[-1, -1]
    solver_in["ubx"]=[1, 0.2]
    solver_out = solver(**solver_in)
    self.assertAlmostEqual(solver_out["x"][0],1,8,"X_opt")
    self.assertAlmostEqual(solver_out["x"][1],0.2,8,"X_opt")

    self.assertAlmostEqual(fmax(solver_out["lam_x"],0)[0],1,8,"Cost should be linear in y0")
    self.assertAlmostEqual(fmax(solver_out["lam_x"],0)[1],(sqrt(p0)*(te*yc0**2-yc0+p0*te)*tan(arctan(yc0/sqrt(p0))+sqrt(p0)*te)+yc0**2)/(2*p0*yc0**2+2*p0**2),8,"Cost should be linear in y0")
    self.assertAlmostEqual(-solver_out["f"][0],(2*y0-log(yc0**2/p0+1))/2-log(cos(arctan(yc0/sqrt(p0))+sqrt(p0)*te)),7,"Cost")
    self.assertAlmostEqual(fmax(-solver_out["lam_x"],0)[0],0,8,"Constraint is supposed to be unactive")
    self.assertAlmostEqual(fmax(-solver_out["lam_x"],0)[1],0,8,"Constraint is supposed to be unactive")

  @requires_nlpsol("ipopt")
  def test_singleshooting2(self):
    self.message("Single shooting 2")
    p0 = 0.2
    y0= 0.2
    yc0=dy0=0.1
    te=0.4

    t=SX.sym("t")
    q=SX.sym("y",2,1)
    p=SX.sym("p",1,1)
    # y
    # y'
    dae={'x':q, 'p':p, 't':t, 'ode':vertcat(*[q[1],p[0]+q[1]**2 ])}
    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    opts["verbose"] = False
    opts["steps_per_checkpoint"] = 10000
    opts["t0"] = 0
    opts["tf"] = te
    integrator = casadi.integrator("integrator", "cvodes", dae, opts)

    var = MX.sym("var",2,1)
    par = MX.sym("par",1,1)

    q0   = vertcat(*[var[0],par])
    parl  = var[1]
    qend = integrator(x0=q0,p=parl)["xf"]

    parc = MX(dy0)

    f = Function('f', [var,par],[qend[0]])
    nlp = {'x':var, 'f':-f(var,parc), 'g':var[0]-var[1]}
    opts = {}
    opts["ipopt.tol"] = 1e-12
    opts["ipopt.hessian_approximation"] = "limited-memory"
    opts["ipopt.max_iter"] = 10
    opts["ipopt.derivative_test"] = "first-order"
    #opts["ipopt.print_level"] = 0
    solver = nlpsol("solver", "ipopt", nlp, opts)
    solver_in = {}
    solver_in["lbx"]=[-1, -1]
    solver_in["ubx"]=[1, 0.2]
    solver_in["lbg"]=[-1]
    solver_in["ubg"]=[0]
    solver_out = solver(**solver_in)

    self.assertAlmostEqual(solver_out["x"][0],0.2,6,"X_opt")
    self.assertAlmostEqual(solver_out["x"][1],0.2,6,"X_opt")

    self.assertAlmostEqual(fmax(solver_out["lam_x"],0)[0],0,8,"Constraint is supposed to be unactive")
    dfdp0 = (sqrt(p0)*(te*yc0**2-yc0+p0*te)*tan(arctan(yc0/sqrt(p0))+sqrt(p0)*te)+yc0**2)/(2*p0*yc0**2+2*p0**2)
    self.assertAlmostEqual(fmax(solver_out["lam_x"],0)[1],1+dfdp0,8)
    self.assertAlmostEqual(solver_out["lam_g"][0],1,8)
    self.assertAlmostEqual(-solver_out["f"][0],(2*y0-log(yc0**2/p0+1))/2-log(cos(arctan(yc0/sqrt(p0))+sqrt(p0)*te)),7,"Cost")
    self.assertAlmostEqual(fmax(-solver_out["lam_x"],0)[0],0,8,"Constraint is supposed to be unactive")
    self.assertAlmostEqual(fmax(-solver_out["lam_x"],0)[1],0,8,"Constraint is supposed to be unactive")

  @requiresPlugin(XmlFile,"tinyxml")
  def test_XML(self):
    self.message("JModelica XML parsing")
    ivp = DaeBuilder()
    ivp.parse_fmi('data/cstr.xml')

    # Separate differential and algebraic variables
    ivp.split_dae()

    self.assertTrue(len(ivp.q)==0)
    self.assertTrue(len(ivp.y)==1)
    m = vertcat(*ivp.ydef)
    self.assertTrue(isinstance(m,MX))
    self.assertEqual(str(m),'cost')
    print(dir(ivp))
    self.assertEqual(len(ivp.dae),3)
    print(type(ivp.s))
    self.assertEqual(len(ivp.s),3) # there are three states
    c = ivp("cstr.c")
    T = ivp("cstr.T")
    cost = ivp("cost")
    self.assertTrue(isinstance(c,MX))

    self.assertEqual(c.name(),"cstr.c")
    self.assertEqual(T.name(),"cstr.T")
    self.assertEqual(cost.name(),"cost")
    self.assertEqual(ivp.nominal("cstr.c"),1000)

    u = ivp("u")
    #self.assertEquals(ivp.path.nnz(),3)
    #self.assertEquals(len(ivp.cfcn_lb),3)
    #self.assertEquals(len(ivp.cfcn_ub),3)
    #self.assertTrue(ivp.cfcn[0].is_equal(T))
    #self.assertTrue(ivp.cfcn[1].is_equal(u))
    #self.assertTrue(ivp.cfcn[2].is_equal(u))
    #self.assertTrue(ivp.cfcn_lb[0].isMinusInf())
    #self.assertEquals(ivp.cfcn_lb[1].to_double(),230)
    #self.assertTrue(ivp.cfcn_lb[2].isMinusInf())
    #self.assertEquals(ivp.cfcn_ub[0].to_double(),350)
    #self.assertTrue(ivp.cfcn_ub[1].isInf())
    #self.assertEquals(ivp.cfcn_ub[2].to_double(),370)
    print(ivp.init)
    print(c,T,cost)
    #print c.atTime(0)
    f=Function('f', [vertcat(*[c,T,cost])],[vertcat(*ivp.init)])
    return
    f_out = f(f_in)
    self.checkarray(f_out[0],matrix([-956.271065,-250.051971,0]).T,"initeq")


    mystates = []

if __name__ == '__main__':
    unittest.main()

