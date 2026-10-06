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
from numpy import inf, pi
from typing import Any, Dict
import casadi as c
import numpy
import unittest
from types import *
from helpers import args, casadiTestCase, codegen_check_digits, memory_heavy, requires_nlpsol

import os

# ipmc codegen links against IPMC_ROOT/build/ipmc and needs casadi's include tree
IPMC_CODEGEN = False
if "SKIP_IPMC_TESTS" not in os.environ and ca.has_nlpsol("ipmc"):
  _ipmc_root = os.environ.get("IPMC_ROOT", "/missing")
  _ipmc_libdir = os.path.join(_ipmc_root, "build", "ipmc")
  if not os.path.isdir(_ipmc_libdir):
    print("IPMC_ROOT not set or has no build/ipmc, skipping ipmc codegen checks")
  elif not os.path.isdir(ca.GlobalOptions.getCasadiIncludePath()):
    print("no casadi include tree at %s, skipping ipmc codegen checks"
          % ca.GlobalOptions.getCasadiIncludePath())
  else:
    IPMC_CODEGEN = {"std": "c99", "extralibs": ["ipmc", "blasfeo"],
                    "extralibdirs": [_ipmc_libdir], "extra_include": [_ipmc_root],
                    "extra_options": [] if os.name == 'nt' else ["-Wno-strict-prototypes"]}

class OCPtests(casadiTestCase):

  def fatrop_case(self,N=2, nx0=2, nu0=2, nx1=2, nu1=2, nx2=2, nu2=2, ng1=2, ng2=2, ng3=2, sp=None,eq=None):
        print("fatrop_case",N,nx0,nu0,nx1,nu1,nx2,nu2,ng1,ng2,ng3,sp,eq)
        if sp is None:
            sp = {}
        if eq is None:
            eq = set()
        nx = [nx0 ,nx1, nx2]
        nu = [nu0, nu1, nu2]
        ng = [ng1, ng2, ng3]
        
        print("nx",nx)
        print("nu",nu)
        print("ng",ng)
        
        ca.DM.rng(1)
        
        A0 = ca.DM.rand(nx1, nx0)
        if "A0" in sp: A0 = ca.project(A0, sp["A0"])
        B0 = ca.DM.rand(nx1, nu0)
        if "B0" in sp: B0 = ca.project(B0, sp["B0"])
        C0 = ca.DM.rand(ng1, nx0)
        if "C0" in sp: C0 = ca.project(C0, sp["C0"])
        D0 = ca.DM.rand(ng1, nu0)
        I0 = ca.DM.eye(nx1)
        
        A1 = ca.DM.rand(nx2, nx1)
        if "A1" in sp: A1 = ca.project(A1, sp["A1"])
        B1 = ca.DM.rand(nx2, nu1)
        if "B1" in sp: B1 = ca.project(B1, sp["B1"])
        C1 = ca.DM.rand(ng2, nx1)
        if "C1" in sp: C1 = ca.project(C1, sp["C1"])
        D1 = ca.DM.rand(ng2, nu1)
        if "D1" in sp: D1 = ca.project(D1, sp["D1"])
        I1 = ca.DM.eye(nx2)
       
        C2 = ca.DM.rand(ng3, nx2)
        if "C2" in sp: C2 = ca.project(C2, sp["C2"])
        D2 = ca.DM.rand(ng3, nu2)
        if "D2" in sp: D2 = ca.project(D2, sp["D2"])
        
        A = ca.blockcat([[A0,B0,I0,ca.DM(nx1,nu1+nx2+nu2)],[C0,D0,ca.DM(ng1,nx1+nu1+nx2+nu2)],[ca.DM(nx2,nx0+nu0),A1,B1,I1,ca.DM(nx2,nu2)],[ca.DM(ng2,nx0+nu0),C1,D1,ca.DM(ng2,nx2+nu2)],[ca.DM(ng3,nx0+nu0+nx1+nu1),C2,D2]])
        
        
        
        equality = [True]*nx1+["ng1" in eq]*ng1+[True]*nx2+["ng2" in eq]*ng2+["ng3" in eq]*ng3
        
        A.sparsity().spy()
        print(A)
       
        x0 = ca.MX.sym("x0",nx0)
        u0 = ca.MX.sym("u0",nu0)
        x1 = ca.MX.sym("x1",nx1)
        u1 = ca.MX.sym("u1",nu1)
        x2 = ca.MX.sym("x2",nx2)
        u2 = ca.MX.sym("u2",nu2)      
        
        x = ca.vertcat(x0,u0,x1,u1,x2,u2)
        nlp = {}
        nlp["x"] = x
        nlp["g"] = ca.DM.zeros(A.shape[0],1) + A @ x
        
        nlp["f"] = ca.sumsqr(x-ca.DM.rand(x.numel(),1))
        
        a = 10
        lbg = ca.vertcat(ca.DM.zeros(nx1,1),-a*ca.DM.ones(ng1,1),ca.DM.zeros(nx2,1),-a*ca.DM.ones(ng2,1),-a*ca.DM.ones(ng3,1))
        ubg = ca.vertcat(ca.DM.zeros(nx1,1),a*ca.DM.ones(ng1,1),ca.DM.zeros(nx2,1),a*ca.DM.ones(ng2,1),a*ca.DM.ones(ng3,1))
        
        print(lbg)

        
        options = {"structure_detection": "manual", "N":N, "nx": nx, "nu":nu, "ng": ng, "equality": equality,"fatrop":{"tol":1e-7}}
        solver = ca.nlpsol("solver","fatrop",nlp,options)
        sol = solver(lbg=lbg,ubg=ubg)

        solver = ca.nlpsol("solver","fatrop",nlp,{"structure_detection": "none", "error_on_fail":True, "equality": equality,"fatrop":{"tol":1e-7}})
        ref = solver(lbg=lbg,ubg=ubg)
        
        for k in sol.keys():
            self.checkarray(sol[k],ref[k],failmessage=k+str(options),digits=6)

        options = {"structure_detection": "auto", "debug":True, "equality": equality,"fatrop":{"tol":1e-7}}
        print(options)
        solver = ca.nlpsol("solver","fatrop",nlp,options)
        sol = solver(lbg=lbg,ubg=ubg)
        
        stats = solver.stats()
        print(stats)
        if nx2>0:
            self.assertTrue(stats["N"]>1)
        
        
        for k in sol.keys():
            self.checkarray(sol[k],ref[k],failmessage=k+str(options),digits=6)

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

    X=ca.SX.sym("X",N+1)
    U=ca.SX.sym("U",N)

    V = ca.vertcat(*[X,U])

    cost = 0
    for i in range(N):
      cost = cost + s*X[i]**2+r*U[i]**2
    cost = cost + q*X[N]**2

    nlp = {'x':V, 'f':cost, 'g':ca.vertcat(*[X[0]-x0,X[1:,0]-(a*X[:N,0]+b*U)])}
    opts = {}
    opts["ipopt.tol"] = 1e-5
    opts["ipopt.hessian_approximation"] = "limited-memory"
    opts["ipopt.max_iter"] = 100
    opts["ipopt.print_level"] = 0
    solver = ca.nlpsol("solver", "ipopt", nlp, opts)
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

    t=ca.SX.sym("t")
    q=ca.SX.sym("y",2,1)
    p=ca.SX.sym("p",1,1)
    # y
    # y'
    dae={'x':q, 'p':p, 't':t, 'ode':ca.vertcat(*[q[1],p[0]+q[1]**2 ])}
    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    opts["verbose"] = False
    opts["steps_per_checkpoint"] = 10000
    integrator = ca.integrator("integrator", "cvodes", dae, 0, te, opts)

    var = ca.MX.sym("var",2,1)
    par = ca.MX.sym("par",1,1)
    parMX= par

    q0   = ca.vertcat(*[var[0],par])
    par  = var[1]
    qend = integrator(x0=q0, p=par)["xf"]

    parc = ca.MX(0)

    f = ca.Function('f', [var,parMX],[qend[0]])
    nlp = {'x':var, 'f':-f(var,parc)}
    opts = {}
    opts["ipopt.tol"] = 1e-12
    opts["ipopt.hessian_approximation"] = "limited-memory"
    opts["ipopt.max_iter"] = 10
    opts["ipopt.derivative_test"] = "first-order"
    opts["ipopt.print_level"] = 0
    solver = ca.nlpsol("solver", "ipopt", nlp, opts)
    solver_in = {}
    solver_in["lbx"]=[-1, -1]
    solver_in["ubx"]=[1, 0.2]
    solver_out = solver(**solver_in)
    print(solver_out["x"])
    self.assertAlmostEqual(solver_out["x"][0],1,7,"X_opt")
    self.assertAlmostEqual(solver_out["x"][1],0.2,7,"X_opt")
    self.assertAlmostEqual(ca.fmax(solver_out["lam_x"],0)[0],1,8,"Cost should be linear in y0")
    self.assertAlmostEqual(ca.fmax(solver_out["lam_x"],0)[1],(ca.sqrt(p0)*(te*yc0**2-yc0+p0*te)*ca.tan(ca.arctan(yc0/ca.sqrt(p0))+ca.sqrt(p0)*te)+yc0**2)/(2*p0*yc0**2+2*p0**2),8,"Cost should be linear in y0")
    self.assertAlmostEqual(-solver_out["f"][0],(2*y0-ca.log(yc0**2/p0+1))/2-ca.log(ca.cos(ca.arctan(yc0/ca.sqrt(p0))+ca.sqrt(p0)*te)),7,"Cost")
    self.assertAlmostEqual(ca.fmax(-solver_out["lam_x"],0)[0],0,8,"Constraint is supposed to be unactive")
    self.assertAlmostEqual(ca.fmax(-solver_out["lam_x"],0)[1],0,8,"Constraint is supposed to be unactive")

  @requires_nlpsol("ipopt")
  def test_singleshooting2(self):
    self.message("Single shooting 2")
    p0 = 0.2
    y0= 0.2
    yc0=dy0=0.1
    te=0.4

    t=ca.SX.sym("t")
    q=ca.SX.sym("y",2,1)
    p=ca.SX.sym("p",1,1)
    # y
    # y'
    dae={'x':q, 'p':p, 't':t, 'ode':ca.vertcat(*[q[1],p[0]+q[1]**2 ])}
    opts = {}
    opts["reltol"] = 1e-15
    opts["abstol"] = 1e-15
    opts["verbose"] = False
    opts["steps_per_checkpoint"] = 10000
    integrator = ca.integrator("integrator", "cvodes", dae, 0, te, opts)

    var = ca.MX.sym("var",2,1)
    par = ca.MX.sym("par",1,1)

    q0   = ca.vertcat(*[var[0],par])
    parl  = var[1]
    qend = integrator(x0=q0,p=parl)["xf"]

    parc = ca.MX(dy0)

    f = ca.Function('f', [var,par],[qend[0]])
    nlp = {'x':var, 'f':-f(var,parc), 'g':var[0]-var[1]}
    opts = {}
    opts["ipopt.tol"] = 1e-12
    opts["ipopt.hessian_approximation"] = "limited-memory"
    opts["ipopt.max_iter"] = 10
    opts["ipopt.derivative_test"] = "first-order"
    #opts["ipopt.print_level"] = 0
    solver = ca.nlpsol("solver", "ipopt", nlp, opts)
    solver_in = {}
    solver_in["lbx"]=[-1, -1]
    solver_in["ubx"]=[1, 0.2]
    solver_in["lbg"]=[-1]
    solver_in["ubg"]=[0]
    solver_out = solver(**solver_in)

    self.assertAlmostEqual(solver_out["x"][0],0.2,6,"X_opt")
    self.assertAlmostEqual(solver_out["x"][1],0.2,6,"X_opt")

    self.assertAlmostEqual(ca.fmax(solver_out["lam_x"],0)[0],0,8,"Constraint is supposed to be unactive")
    dfdp0 = (ca.sqrt(p0)*(te*yc0**2-yc0+p0*te)*ca.tan(ca.arctan(yc0/ca.sqrt(p0))+ca.sqrt(p0)*te)+yc0**2)/(2*p0*yc0**2+2*p0**2)
    self.assertAlmostEqual(ca.fmax(solver_out["lam_x"],0)[1],1+dfdp0,8)
    self.assertAlmostEqual(solver_out["lam_g"][0],1,8)
    self.assertAlmostEqual(-solver_out["f"][0],(2*y0-ca.log(yc0**2/p0+1))/2-ca.log(ca.cos(ca.arctan(yc0/ca.sqrt(p0))+ca.sqrt(p0)*te)),7,"Cost")
    self.assertAlmostEqual(ca.fmax(-solver_out["lam_x"],0)[0],0,8,"Constraint is supposed to be unactive")
    self.assertAlmostEqual(ca.fmax(-solver_out["lam_x"],0)[1],0,8,"Constraint is supposed to be unactive")

  @requires_nlpsol("fatrop")
  @requires_nlpsol("ipopt")
  @memory_heavy()
  def test_fatrop(self):
  
  
    flags = []
    if os.name != 'nt':
      flags = ["-Wno-strict-prototypes"]
  
    def test_problems():
    
        for i in range(2):

            T = 10. # Time horizon
            N = 10 # number of control intervals

            # Declare model variables
            x1 = ca.MX.sym('x1')
            x2 = ca.MX.sym('x2')
            x = ca.vertcat(x1, x2)
            u = ca.MX.sym('u')
            p = ca.MX.sym('p')

            # Model equations
            xdot = ca.vertcat((1-x2**2)*x1 - x2 + u+p, x1)

            F = ca.integrator("F","rk",{"x":x,"p":p,"u":u,"ode":xdot}, 0, 1, {"simplify":True,"number_of_finite_elements":1})

            # Start with an empty NLP
            w=[]
            w0 = []
            lbw = []
            ubw = []
            J = 0
            g=[]
            lbg = []
            ubg = []
            equality = []

            # "Lift" initial conditions
            Xk = ca.MX.sym('X0', 2)
            w += [Xk]
            lbw += [0, 1]
            ubw += [0, 1]
            w0 += [0.1, 0.2]

            # Formulate the NLP
            for k in range(N):
                # New NLP variable for the control
                Uk = ca.MX.sym('U_' + str(k))
                w   += [Uk]
                lbw += [-1]
                ubw += [1]
                w0  += [0.3]

                # Integrate till the end of the interval
                Fk = F(x0=Xk, u=Uk, p=p)
                Xk_end = Fk['xf']
                J=J+ca.sumsqr(Xk)+ca.sumsqr(Uk)



                # New NLP variable for state at end of interval
                Xk_next = ca.MX.sym('X_' + str(k+1), 2)
                w   += [Xk_next]
                lbw += [-0.25 if i==0 else -inf, -inf]
                ubw += [  inf,  inf]
                w0  += [0.1, 0.2]
                    
                # Add equality constraint
                g   += [Xk_next-Xk_end]
                lbg += [0, 0]
                ubg += [0, 0]
                equality+= [True,True]

                if i>=1:
                    g   += [ca.sin(Xk[0])]
                    lbg += [-0.25]
                    ubg += [inf]
                    equality+= [False]
                    
                Xk = Xk_next
            if i>=2:
                    
                # "Lift" initial conditions
                Xk = ca.MX.sym('X0', 2)
                w += [Xk]
                lbw += [-inf, -inf]
                ubw += [inf, inf]
                w0 += [0.1, 0.2]
                
                
                # Add equality constraint
                g   += [Xk_next-Xk]
                lbg += [0, 0]
                ubg += [0, 0]
                equality+= [True,True]

                # Formulate the NLP
                for k in range(N):
                    # New NLP variable for the control
                    Uk = ca.MX.sym('U_' + str(k))
                    w   += [Uk]
                    lbw += [-0.1]
                    ubw += [0.1]
                    w0  += [0.3]

                    # Integrate till the end of the interval
                    Fk = F(x0=Xk, u=Uk, p=p)
                    Xk_end = Fk['xf']
                    J=J+3*ca.sumsqr(Xk)+ca.sumsqr(Uk)

                    # New NLP variable for state at end of interval
                    Xk_next = ca.MX.sym('X_' + str(k+1), 2)
                    w   += [Xk_next]
                    lbw += [-inf, -inf]
                    ubw += [  inf,  inf]
                    w0  += [0.1, 0.2]
                        
                    # Add equality constraint
                    g   += [Xk_next-Xk_end]
                    lbg += [0, 0]
                    ubg += [0, 0]
                    equality+= [True,True]

                    Xk = Xk_next
                    
            if i>=3:
                # Declare model variables
                x = ca.MX.sym('x',3)
                u = ca.MX.sym('u',2)
                
                A = ca.DM([[1,0,0.3],[0,1,0.7],[0.2,0,1]])
                B = ca.DM([[1,0],[0,1],[0.5,0.5]])
                D = ca.DM([[0.2,0.3],[0.8,0.7],[0.1,1]])

                F = ca.Function("F",[x,u],[A @ x+B @ u])

                    
                # "Lift" initial conditions
                Xk = ca.MX.sym('X0', 3)
                w += [Xk]
                lbw += [-inf, -inf, -inf]
                ubw += [inf, inf, inf]
                w0 += [0.1, 0.2, 0.3]
                
                
                # Add equality constraint
                g   += [D @ Xk_next-Xk]
                lbg += [0, 0, 0]
                ubg += [0, 0, 0]
                equality+= [True,True, True]

                # Formulate the NLP
                for k in range(N):
                    # New NLP variable for the control
                    Uk = ca.MX.sym('U_' + str(k),2)
                    w   += [Uk]
                    lbw += [-1,-1]
                    ubw += [1,1]
                    w0  += [0.3,0.3]

                    # Integrate till the end of the interval
                    Xk_end = F(Xk, Uk)
                    J=J+ca.sumsqr(Xk)+ca.sumsqr(Uk)

                    # New NLP variable for state at end of interval
                    Xk_next = ca.MX.sym('X_' + str(k+1), 3)
                    w   += [Xk_next]
                    lbw += [-inf, -inf, -inf]
                    ubw += [  inf,  inf, inf]
                    w0  += [0.1, 0.2, 0.3]
                        
                    # Add equality constraint
                    g   += [Xk_next-Xk_end]
                    lbg += [0, 0, 0]
                    ubg += [0, 0, 0]
                    equality+= [True,True, True]

                    Xk = Xk_next
 
                        
                
            # Create an NLP solver
            yield {'f': J, 'x': ca.vertcat(*w), 'g': ca.vertcat(*g), 'p': p}, dict(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg, p=0),equality
        
        for i in range(1):
            # Multi-stage with varying number of inequalities

            T = 10. # Time horizon
            N = 10 # number of control intervals

            # Declare model variables
            x1 = ca.MX.sym('x1')
            x2 = ca.MX.sym('x2')
            x = ca.vertcat(x1, x2)
            u = ca.MX.sym('u')
            p = ca.MX.sym('p')

            # Model equations
            xdot = ca.vertcat((1-x2**2)*x1 - x2 + u+p, x1)

            F = ca.integrator("F","rk",{"x":x,"p":p,"u":u,"ode":xdot}, 0, 1, {"simplify":True,"number_of_finite_elements":1})

            # Start with an empty NLP
            w=[]
            w0 = []
            lbw = []
            ubw = []
            J = 0
            g=[]
            lbg = []
            ubg = []
            equality = []

            # "Lift" initial conditions
            Xk = ca.MX.sym('X0', 2)
            w += [Xk]
            lbw += [0, 1]
            ubw += [0, 1]
            w0 += [0.1, 0.2]

            # Formulate the NLP
            for k in range(N):
                # New NLP variable for the control
                Uk = ca.MX.sym('U1_' + str(k))
                w   += [Uk]
                lbw += [-1]
                ubw += [1]
                w0  += [0.3]

                # Integrate till the end of the interval
                Fk = F(x0=Xk, u=Uk, p=p)
                Xk_end = Fk['xf']
                J=J+ca.sumsqr(Xk)+ca.sumsqr(Uk)

                # New NLP variable for state at end of interval
                Xk_next = ca.MX.sym('X1_' + str(k+1), 2)
                w   += [Xk_next]
                lbw += [-0.25, -inf]
                ubw += [  inf,  inf]
                w0  += [0.1, 0.2]
                    
                # Add equality constraint
                g   += [Xk_next-Xk_end]
                lbg += [0, 0]
                ubg += [0, 0]
                equality += [True, True]

                Xk = Xk_next

            # New NLP variable for the control
            Uk = ca.MX.sym('U1')
            w   += [Uk]
            lbw += [-1]
            ubw += [1]
            w0  += [0.3]

            # Integrate till the end of the interval
            Fk = F(x0=Xk, u=Uk, p=p)
            Xk_end = Fk['xf']
            J=J+Xk[0]**2


            # New NLP variable for state at end of interval
            Xk = ca.MX.sym('X1', 2)
            w   += [Xk]
            lbw += [-inf, -inf]
            ubw += [  inf,  inf]
            w0  += [0.1, 0.2]
                
            # Add equality constraint
            g   += [Xk-Xk_end]
            lbg += [0, 0]
            ubg += [0, 0]
            equality += [True, True]


            # "Lift" initial conditions
            #Xk = MX.sym('X0', 2)
            #w += [Xk]
            #lbw += [-inf, -inf]
            #ubw += [inf, inf]
            #w0 += [0.1, 0.2]


            # Add equality constraint
            #g   += [Xk_next-Xk]
            #lbg += [0, 0]
            #ubg += [0, 0]

            A = ca.DM([[1,0.1],[0.2,1.1]])
            B = ca.DM([[0.2],[0.7]])

            F = ca.Function("F",[x,u],[A @ x+B @ u])

            # Formulate the NLP
            for k in range(N):
                # New NLP variable for the control
                Uk = ca.MX.sym('U2_' + str(k))
                w   += [Uk]
                lbw += [-inf]
                ubw += [inf]
                w0  += [0.3]


                # Integrate till the end of the interval
                Xk_end = F(Xk, Uk)
                J=J+3*ca.sumsqr(Xk)+ca.sumsqr(Uk)

                # New NLP variable for state at end of interval
                Xk_next = ca.MX.sym('X2_' + str(k+1), 2)
                w   += [Xk_next]
                lbw += [-inf, -inf]
                ubw += [  inf,  inf]
                w0  += [0.1, 0.2]
                    
                # Add equality constraint
                g   += [Xk_next-Xk_end]
                lbg += [0, 0]
                ubg += [0, 0]
                equality += [True, True]

                g   += [2*Uk]
                lbg += [-0.1]
                ubg += [0.1]
                equality += [False]


                Xk = Xk_next
        
        
            # Create an NLP solver
            yield {'f': J, 'x': ca.vertcat(*w), 'g': ca.vertcat(*g), 'p': p}, dict(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg, p=0), equality
        
        for i in range(1):
            # Equality constraints

            T = 10. # Time horizon
            N = 10 # number of control intervals

            # Declare model variables
            x1 = ca.MX.sym('x1')
            x2 = ca.MX.sym('x2')
            x = ca.vertcat(x1, x2)
            u1 = ca.MX.sym('u1')
            u2 = ca.MX.sym('u2')
            u = ca.vertcat(u1, u2)
            p = ca.MX.sym('p')

            # Model equations
            xdot = ca.vertcat((1-x2**2)*x1 - x2 + u1+p, x1+u2)

            F = ca.integrator("F","rk",{"x":x,"p":p,"u":u,"ode":xdot}, 0, 1, {"simplify":True,"number_of_finite_elements":1})

            # Start with an empty NLP
            w=[]
            w0 = []
            lbw = []
            ubw = []
            J = 0
            g=[]
            lbg = []
            ubg = []
            equality = []

            # "Lift" initial conditions
            Xk = ca.MX.sym('X0', 2)
            w += [Xk]
            lbw += [0, 1]
            ubw += [0, 1]
            w0 += [0.1, 0.2]

            # Formulate the NLP
            for k in range(N):
                # New NLP variable for the control
                Uk = ca.MX.sym('U_' + str(k),2)
                w   += [Uk]
                lbw += [-1,0.1]
                ubw += [1,0.1]
                w0  += [0.3,0]

                # Integrate till the end of the interval
                Fk = F(x0=Xk, u=Uk, p=p)
                Xk_end = Fk['xf']
                J=J+ca.sumsqr(Xk)+ca.sumsqr(Uk)



                # New NLP variable for state at end of interval
                Xk_next = ca.MX.sym('X_' + str(k+1), 2)
                w   += [Xk_next]
                lbw += [-0.25 if i==0 else -inf, -inf]
                ubw += [  inf,  inf]
                w0  += [0.1, 0.2]
                    
                # Add equality constraint
                g   += [Xk_next-Xk_end]
                lbg += [0, 0]
                ubg += [0, 0]
                equality += [True,True]

                Xk = Xk_next
            # Create an NLP solver
            yield {'f': J, 'x': ca.vertcat(*w), 'g': ca.vertcat(*g), 'p': p}, dict(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg, p=0), equality
            

        T = 10. # Time horizon
        N = 10 # number of control intervals

        # Declare model variables
        x1 = ca.MX.sym('x1')
        x2 = ca.MX.sym('x2')
        x = ca.vertcat(x1, x2)
        u = ca.MX.sym('u')
        p = ca.MX.sym('p')

        # Model equations
        xdot = ca.vertcat((1-x2**2)*x1 - x2 + u+p, x1)

        F = ca.integrator("F","rk",{"x":x,"p":p,"u":u,"ode":xdot}, 0, 1, {"simplify":True,"number_of_finite_elements":1})

        # Start with an empty NLP
        w=[]
        w0 = []
        lbw = []
        ubw = []
        J = 0
        g=[]
        lbg = []
        ubg = []
        equality = []

        # "Lift" initial conditions
        Xk = ca.MX.sym('X0', 2)
        w += [Xk]
        lbw += [0, 1]
        ubw += [0, 1]
        w0 += [0.1, 0.2]

        # Formulate the NLP
        for k in range(N):
            # New NLP variable for the control
            Uk = ca.MX.sym('U1_' + str(k))
            w   += [Uk]
            lbw += [-1]
            ubw += [1]
            w0  += [0.3]

            # Integrate till the end of the interval
            Fk = F(x0=Xk, u=Uk, p=p)
            Xk_end = Fk['xf']
            J=J+ca.sumsqr(Xk)+ca.sumsqr(Uk)

            # New NLP variable for state at end of interval
            Xk_next = ca.MX.sym('X1_' + str(k+1), 2)
            w   += [Xk_next]
            lbw += [-0.25, -inf]
            ubw += [  inf,  inf]
            w0  += [0.1, 0.2]
                
            # Add equality constraint
            g   += [Xk_next-Xk_end]
            lbg += [0, 0]
            ubg += [0, 0]
            equality += [True,True]

            Xk = Xk_next

        J=J+Xk[0]**2

        A = ca.DM([[1,0,0.3],[0,1,0.7],[0.2,0,1]])
        B = ca.DM([[1,0],[0,1],[0.5,0.5]])
        D = ca.DM([[0.2,0.3],[0.8,0.7],[0.1,1]])

        # New NLP variable for state at end of interval
        Xk = ca.MX.sym('X1', 3)
        w   += [Xk]
        lbw += [-inf, -inf, -inf]
        ubw += [  inf,  inf, inf]
        w0  += [0.7, 0.8, 0.9]
            
        # Add equality constraint
        g   += [Xk-D @ Xk_next]
        lbg += [0, 0, 0]
        ubg += [0, 0, 0]
        equality += [True,True,True]

        u = ca.MX.sym("u",2)
        x = ca.MX.sym("x",3)
        F = ca.Function("F",[x,u],[A @ x+B @ u])

        # Formulate the NLP
        for k in range(N):
            # New NLP variable for the control
            Uk = ca.MX.sym('U2_' + str(k),2)
            w   += [Uk]
            lbw += [-0.1,-0.1]
            ubw += [0.1,0.1]
            w0  += [1.3,1.2]


            # Integrate till the end of the interval
            Xk_end = F(Xk, Uk)
            J=J+3*ca.sumsqr(Xk)+ca.sumsqr(Uk)

            # New NLP variable for state at end of interval
            Xk_next = ca.MX.sym('X2_' + str(k+1), 3)
            w   += [Xk_next]
            lbw += [-inf, -inf, -inf]
            ubw += [  inf,  inf, inf]
            w0  += [0.7, 0.8, 0.9]
                
            # Add equality constraint
            g   += [Xk_next-Xk_end]
            lbg += [0, 0, 0]
            ubg += [0, 0, 0]
            equality += [True,True,True]


            Xk = Xk_next
 
         # Create an NLP solver
        yield {'f': J, 'x': ca.vertcat(*w), 'g': ca.vertcat(*g), 'p': p}, dict(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg, p=0), equality
        
    local_codegen_check_digits = codegen_check_digits
    if os.name == 'nt': local_codegen_check_digits = local_codegen_check_digits -1
    for i,(prob,args,equality) in enumerate(test_problems()):
    
        ca.jacobian_sparsity(prob["g"],prob["x"]).spy()

        solutions = {}
        stats = {}
        # fixed_variable_treatment: the horizon pins x_0 with lbw == ubw, and
        # ipopt's default ("make_parameter") drops such variables from the NLP
        # and reports lam_x = 0 for them, while fatrop reports the real
        # multiplier (-2.7, -10.2 on the first problem here).  The lam_x
        # comparison below is only meaningful if ipopt is asked to keep them.
        for solver, solver_options in [("ipopt",{"ipopt":{"fixed_variable_treatment":"make_constraint"}}),("fatrop",{"structure_detection":"auto","fatrop":{"tol":1e-8,"max_iter":100},"equality":equality})]:
            f = ca.nlpsol('solver', solver, prob, solver_options)
            #if solver=="fatrop" and i==2: raise Exception() 

            # Solve the NLP
            solutions[solver] = f(**args)
            stats[solver] = f.stats()
            
            if solver!="ipopt":
                self.check_codegen(f,args,std="c99",extralibs=["fatrop","blasfeo"],extra_options=flags,digits=local_codegen_check_digits)
                self.check_serialize(f,args)
        
        for k in solutions["ipopt"].keys():
            if k in ["x","f","g","lam_g","lam_x","lam_p"]:
                v_ref = solutions["ipopt"][k]
                v = solutions["fatrop"][k]
                
                self.checkarray(v,v_ref,failmessage=k,digits=5)
        assert(abs(stats["ipopt"]["iter_count"]-stats["fatrop"]["iter_count"])<=2)



  @requires_nlpsol("fatrop")
  def test_fatrop_sanitize(self):
  
  
    def test_problems():
    

            
            T = 10. # Time horizon
            N = 10 # number of control intervals

            # Declare model variables
            x1 = ca.MX.sym('x1')
            x2 = ca.MX.sym('x2')
            x = ca.vertcat(x1, x2)
            u = ca.MX.sym('u')
            p = ca.MX.sym('p')

            # Model equations
            xdot = ca.vertcat((1-x2**2)*x1 - x2 + u+p, x1)

            F = ca.integrator("F","rk",{"x":x,"p":p,"u":u,"ode":xdot}, 0, 1, {"simplify":True,"number_of_finite_elements":1})

            # Start with an empty NLP
            w=[]
            w0 = []
            lbw = []
            ubw = []
            J = 0
            g=[]
            equality = []
            lbg = []
            ubg = []

            # "Lift" initial conditions
            Xk = ca.MX.sym('X0', 2)
            w += [Xk]
            lbw += [0, 1]
            ubw += [0, 1]
            w0 += [0.1, 0.2]

            # Formulate the NLP
            for k in range(N):
                # New NLP variable for the control
                Uk = ca.MX.sym('U_' + str(k))
                w   += [Uk]
                lbw += [-1]
                ubw += [1]
                w0  += [0.3]

                # Integrate till the end of the interval
                Fk = F(x0=Xk, u=Uk, p=p)
                Xk_end = Fk['xf']
                J=J+ca.sumsqr(Xk)+ca.sumsqr(Uk)



                # New NLP variable for state at end of interval
                Xk_next = ca.MX.sym('X_' + str(k+1), 2)
                w   += [Xk_next]
                lbw += [-0.25, -inf]
                ubw += [  inf,  inf]
                w0  += [0.1, 0.2]
                    
                # Add equality constraint
                g   += [Xk_end-Xk_next]
                lbg += [0, 0]
                ubg += [0, 0]
                equality+= [True, True]

                    
                Xk = Xk_next
           
 
             # Create an NLP solver
            yield {'f': J, 'x': ca.vertcat(*w), 'g': ca.vertcat(*g), 'p': p}, dict(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg, p=0),equality
        
    for i,(prob,args,equality) in enumerate(test_problems()):
    
        ca.jacobian_sparsity(prob["g"],prob["x"]).spy()


        for solver, solver_options in [("fatrop",{"structure_detection": "auto", "verbose": True, "equality": equality})]:
            f = ca.nlpsol('solver', solver, prob, solver_options)
            #if solver=="fatrop" and i==2: raise Exception() 

            # Solve the NLP
            with self.assertInAnyOutput("gap-closing"):
                f(**args)
       
  @requires_nlpsol("fatrop")
  def test_detect(self):

    for nx1 in [2]:
        for nx2 in [2,0]:
            for nu1 in [2,0]:
                for nu2 in [2,0]:
                    for ng1 in [2,0]:
                        for ng2 in [2,0]:
                            for ng3 in [2,0]:
                                print("test_detect",nx1,nx2,nu1,nu2,ng1,ng2,ng3)
                                self.fatrop_case(N=2,nx0=2,nu0=2,nx1=nx1,nu1=nu1,nx2=nx2,nu2=nu2,ng1=ng1,ng2=ng2,ng3=ng3)

  @requires_nlpsol("fatrop")
  def test_detect_adversarial(self):
    
    D2 = ca.sparsify(ca.blockcat([[1,0,0],[1,1,1]])).sparsity()
    self.fatrop_case(nu2=3,sp={"D2": D2})
    
    D2 = ca.sparsify(ca.blockcat([[1,0,0],[1,1,1]])).sparsity()
    C2 = ca.sparsify(ca.blockcat([[0,1],[0,0]])).sparsity()
    self.fatrop_case(nu2=3,sp={"D2": D2, "C2": C2})
    with self.assertInAnyOutput("gap-closing constraints must be like"):
        self.fatrop_case(nu2=3,sp={"D2": D2, "C2": C2},eq={'ng3'}) # Why is this not trig
    
    self.fatrop_case(nu0=0,nx0=1)
    
    self.fatrop_case(nx2=0)
    
    
    #with self.assertInAnyOutput("Gap-closing constraint must depend on a state"):
    self.fatrop_case(nu2=3,sp={"A1": ca.Sparsity(2,2), "B1": ca.Sparsity(2,2)})
    #with self.assertInAnyOutput("Gap-closing constraint must depend on a state"):
    self.fatrop_case(nu2=3,sp={"A1": ca.Sparsity(2,2)})
        
    self.fatrop_case(nx2=1,ng1=0,nu1=0)
    
  @requires_nlpsol("fatrop")
  def test_bug(self):

    x = ca.MX.sym("x")

    for structure_detection in ["none","auto"]:

        opts = {"expand": True, "structure_detection": structure_detection,"equality":[True]}
        
        solver = ca.nlpsol("solver","fatrop",{"x":x,"g":x-1},opts)
        self.assertAlmostEqual(solver(lbg=0,ubg=0)["x"],1,5)

        solver = ca.nlpsol("solver","fatrop",{"x":x,"g":x},opts)
        self.assertAlmostEqual(solver(lbg=1,ubg=1)["x"],1,5)

        solver = ca.nlpsol("solver","fatrop",{"x":x,"g":x-2},opts)
        self.assertAlmostEqual(solver(lbg=3,ubg=3)["x"],5,5)
        
  
  # ipmc on staircase OCPs, against ipopt on the same NLP
  IPMC_OPTS = {"print_level": 0, "tol": 1e-9, "max_iter": 500}

  IPOPT_HARD_REF_OPTS = {"print_time": False,
                         "ipopt": {"print_level": 0, "sb": "yes", "tol": 1e-12,
                                   "constr_viol_tol": 1e-12, "dual_inf_tol": 1e-9,
                                   "compl_inf_tol": 1e-12, "acceptable_tol": 1e-12,
                                   "bound_relax_factor": 0,
                                   "fixed_variable_treatment": "make_constraint",
                                   "max_iter": 3000}}

  def sprint_ocp(self, N=15, v_lo=1.0, v_hi=6.0, p_goal=5.0, T=1.0,
                 u_max=40.0) -> Dict[str, Any]:
    """Minimum-effort sprint OCP, RK4 on p' = v, v' = u - 0.02 v^2, staircase layout"""
    dt = T/N
    xs = ca.SX.sym("xs", 2)
    us = ca.SX.sym("us")
    ode = ca.Function("ode", [xs, us], [ca.vertcat(xs[1], us-0.02*xs[1]*xs[1])])
    k1 = ode(xs, us)
    k2 = ode(xs+dt/2*k1, us)
    k3 = ode(xs+dt/2*k2, us)
    k4 = ode(xs+dt*k3, us)
    F = ca.Function("F", [xs, us], [xs+dt/6*(k1+2*k2+2*k3+k4)])

    X = [ca.SX.sym("x_%d" % k, 2) for k in range(N+1)]
    U = [ca.SX.sym("u_%d" % k) for k in range(N)]

    lbx, ubx = [], []
    for k in range(N+1):
      lbx += [0.0, 0.0] if k == 0 else [-inf, -inf]
      ubx += [0.0, 0.0] if k == 0 else [inf, inf]
      if k < N:
        lbx.append(-u_max)
        ubx.append(u_max)

    var, g, lbg, ubg, ng, equality = [], [], [], [], [], []
    for k in range(N):
      var += [X[k], U[k]]
      g.append(X[k+1]-F(X[k], U[k]))     # gap-closing / dynamics rows
      lbg += [0.0, 0.0]
      ubg += [0.0, 0.0]
      equality += [True, True]
      g.append(X[k][1])                  # the speed band, stage k
      lbg.append(v_lo if k >= 1 else -inf)
      ubg.append(v_hi)
      equality.append(False)
      ng.append(1)
    var.append(X[N])
    g.append(X[N][1])
    lbg.append(v_lo)
    ubg.append(v_hi)
    equality.append(False)
    g.append(X[N][0])                    # terminal reach constraint
    lbg.append(p_goal)
    ubg.append(inf)
    equality.append(False)
    ng.append(2)

    x = ca.vertcat(*var)
    return dict(nlp={"x": x, "f": dt*sum(U[k]*U[k] for k in range(N)),
                     "g": ca.vertcat(*g)},
                bounds=dict(x0=0, lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg),
                equality=equality,
                structure=dict(N=N, nx=[2]*(N+1), nu=[1]*N+[0], ng=ng))

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  def test_ipmc_structured(self):
    self.message("ipmc structure_detection none, manual and auto against ipopt")
    p = self.sprint_ocp()
    ref = ca.nlpsol("reference", "ipopt", p["nlp"],
                    dict(self.IPOPT_HARD_REF_OPTS,
                         ipopt=dict(self.IPOPT_HARD_REF_OPTS["ipopt"])))
    r_ref = ref(**p["bounds"])
    self.assertTrue(ref.stats()["success"])
    # the speed cap binds on most stages, so the duals are compared on active rows
    v = [float(r_ref["x"][3*k+1]) for k in range(1, p["structure"]["N"]+1)]
    n_active = sum(1 for vk in v if vk > 6.0-1e-6)
    print("test_ipmc_structured active speed rows", n_active, "f", float(r_ref["f"]))
    self.assertTrue(n_active >= 5)

    for sd in ["none", "manual", "auto"]:
      print("test_ipmc_structured", sd)
      opts = {"structure_detection": sd, "ipmc": dict(self.IPMC_OPTS)}
      if sd == "manual":
        opts.update(p["structure"])
      elif sd == "auto":
        opts["equality"] = p["equality"]
      solver = ca.nlpsol("solver", "ipmc", p["nlp"], opts)
      r = solver(**p["bounds"])
      self.assertTrue(solver.stats()["success"])

      if sd != "none":
        # the detected horizon is the user's
        st = solver.stats()
        self.assertEqual(st["N"], p["structure"]["N"])
        for k in ["nx", "nu", "ng"]:
          self.checkarray(ca.DM(p["structure"][k]), ca.DM(st[k]),
                          "structure:"+k, digits=12)

      # ipmc relaxes bounds by 1e-8, ipopt does not: digits 4 to 6
      for k, d in [("f", 4), ("x", 5), ("g", 6), ("lam_g", 4), ("lam_x", 4)]:
        self.checkarray(r_ref[k], r[k], sd+":"+k, digits=d)

    # the block partition survives codegen and serialization
    solver = ca.nlpsol("solver", "ipmc", p["nlp"],
                       dict({"structure_detection": "auto",
                             "equality": p["equality"],
                             "ipmc": dict(self.IPMC_OPTS)}))
    if IPMC_CODEGEN:
      self.check_codegen(solver, p["bounds"], **IPMC_CODEGEN)
    self.check_serialize(solver, p["bounds"])

  # nxc: the last nxc states of every stage are constant, x_{k+1} = x_k

  # nxc on vs off differs by at most 4e-14 primal and 2e-11 dual
  NXC_DIGITS = [("f", 12), ("x", 12), ("g", 13), ("lam_g", 10), ("lam_x", 10)]
  # a solve ending in restoration has multipliers of 1e7 to 1e9
  NXC_RESTO_DIGITS = [("f", 6), ("x", 8), ("g", 8), ("lam_g", 6), ("lam_x", 6)]

  def nxc_ocp(self, N=10, nc=4, T=1.0, p_goal=2.0, u_max=40.0, W=10.0,
              bound_every=True, terminal_eq=False) -> Dict[str, Any]:
    """Sprint OCP whose state carries nc constants: in the dynamics, the cost and a path row"""
    dt = T/N
    nx = 2+nc
    xs = ca.SX.sym("xs", nx)
    us = ca.SX.sym("us")
    ode = ca.Function("ode", [xs, us], [ca.vertcat(xs[1], us-xs[2]*xs[1]*xs[1])])
    pad = lambda z: ca.vertcat(z, ca.SX.zeros(nc))
    k1 = ode(xs, us)
    k2 = ode(xs+dt/2*pad(k1), us)
    k3 = ode(xs+dt/2*pad(k2), us)
    k4 = ode(xs+dt*pad(k3), us)
    Fd = ca.Function("Fd", [xs, us], [xs[:2]+dt/6*(k1+2*k2+2*k3+k4)])

    X = [ca.SX.sym("x_%d" % k, nx) for k in range(N+1)]
    U = [ca.SX.sym("u_%d" % k) for k in range(N)]
    Wj = [20.0*(1.0+0.35*j) for j in range(nc)]
    cap = lambda k: X[k][3] + sum([X[k][2+j] for j in range(2, nc)])

    lbx, ubx = [], []
    for k in range(N+1):
      lbx += [0.0, 0.0] if k == 0 else [-inf, -inf]
      ubx += [0.0, 0.0] if k == 0 else [inf, inf]
      # bound_every=False bounds the constants at stage 0 only
      if bound_every or k == 0:
        lbx.append(0.005); ubx.append(0.05)
        if nc > 1:
          lbx.append(0.5); ubx.append(20.0)
        for j in range(2, nc):
          lbx.append(0.0); ubx.append(0.4)
      else:
        lbx += [-inf]*nc; ubx += [inf]*nc
      if k < N:
        lbx.append(-u_max); ubx.append(u_max)

    var, g, lbg, ubg, ng, equality = [], [], [], [], [], []
    # row index of the v_k <= cap row of every stage
    cap_rows = []
    for k in range(N):
      var += [X[k], U[k]]
      g.append(X[k+1][:2]-Fd(X[k], U[k]))            # the real dynamics
      lbg += [0.0]*2; ubg += [0.0]*2; equality += [True]*2
      g.append(X[k+1][2:]-X[k][2:])                  # the [0 I] block
      lbg += [0.0]*nc; ubg += [0.0]*nc; equality += [True]*nc
      if nc > 1:
        cap_rows.append(len(equality))
        g.append(X[k][1]-cap(k))
        lbg.append(-inf); ubg.append(0.0); equality.append(False)
        ng.append(1)
      else:
        ng.append(0)
    var.append(X[N])
    rows = 0
    if terminal_eq:
      # two terminal equalities without an input to pivot on, pushed back to stage N-1
      g.append(X[N][0]); lbg.append(p_goal); ubg.append(p_goal)
      equality.append(True); rows += 1
      if nc > 1:
        g.append(X[N][1]-cap(N)); lbg.append(0.0); ubg.append(0.0)
        equality.append(True); rows += 1
    else:
      if nc > 1:
        cap_rows.append(len(equality))
        g.append(X[N][1]-cap(N)); lbg.append(-inf); ubg.append(0.0)
        equality.append(False); rows += 1
      g.append(X[N][0]); lbg.append(p_goal); ubg.append(inf)
      equality.append(False); rows += 1
    ng.append(rows)

    f = dt*sum(U[k]*U[k] for k in range(N))
    if nc > 1:
      f = f + W*X[0][3]*X[0][3]
    for j in range(2, nc):
      f = f + Wj[j]*X[0][2+j]*X[0][2+j]

    return dict(nlp={"x": ca.vertcat(*var), "f": f, "g": ca.vertcat(*g)},
                bounds=dict(x0=0.1, lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg),
                equality=equality,
                structure=dict(N=N, nx=[nx]*(N+1), nu=[1]*N+[0], ng=ng),
                lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg, cap_rows=cap_rows,
                nxc=nc, nc=nc, N=N)

  def nxc_solve(self, p, sd, nxc=None, opts_extra=None):
    opts = {"structure_detection": sd, "print_time": False, "ipmc": dict(self.IPMC_OPTS)}
    if sd == "manual":
      opts.update(p["structure"])
    elif sd == "auto":
      opts["equality"] = p["equality"]
    if nxc is not None:
      opts["nxc"] = nxc
    if opts_extra:
      opts.update(opts_extra)
    solver = ca.nlpsol("solver", "ipmc", p["nlp"], opts)
    return solver, solver(**p["bounds"])

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  def test_ipmc_nxc(self):
    self.message("ipmc nxc does not change the answer")
    p = self.nxc_ocp(N=10, nc=4)
    N, nc, nx = p["N"], p["nc"], 2+p["nc"]

    ref = ca.nlpsol("reference", "ipopt", p["nlp"],
                    dict(self.IPOPT_HARD_REF_OPTS,
                         ipopt=dict(self.IPOPT_HARD_REF_OPTS["ipopt"])))
    r_ref = ref(**p["bounds"])
    self.assertTrue(ref.stats()["success"])

    # the constants are active at the solution
    xs = list(numpy.array(r_ref["x"]).ravel())
    c = xs[2:2+nc]
    v = [xs[k*(nx+1)+1] for k in range(N+1)]
    capv = c[1] + sum(c[2:])
    n_cap = sum(1 for vk in v if vk > capv - 1e-7)
    print("test_ipmc_nxc constants", [round(ci, 6) for ci in c],
          "active cap rows", n_cap, "f", float(r_ref["f"]))
    # the path constraint that the constants sit in is active on most stages
    self.assertTrue(n_cap >= 5)
    # c_0, the one in the dynamics, is pinned on its lower bound
    self.assertTrue(abs(c[0]-0.005) < 1e-7)
    # one other constant sits on its upper bound and one strictly inside
    self.assertTrue(any(abs(cj-0.4) < 1e-7 for cj in c[2:]))
    self.assertTrue(any(1e-6 < cj < 0.4-1e-6 for cj in c[2:]))
    # the terminal reach is active, which is what pays for the effort
    self.assertTrue(abs(xs[N*(nx+1)] - 2.0) < 1e-7)

    for sd in ["manual", "auto"]:
      s0, r0 = self.nxc_solve(p, sd, None)
      sz, rz = self.nxc_solve(p, sd, 0)
      s1, r1 = self.nxc_solve(p, sd, p["nxc"])
      for s in [s0, sz, s1]:
        self.assertTrue(s.stats()["success"])
      self.assertEqual(s1.stats()["nxc"], p["nxc"])
      self.assertEqual(s0.stats()["nxc"], 0)

      # nxc=0 is bit-identical to nxc absent
      for k in ["f", "x", "g", "lam_g", "lam_x"]:
        self.assertEqual(float(ca.norm_inf(r0[k]-rz[k])), 0.0,
                         sd+": nxc=0 is not bit-identical on "+k)

      # the same iterate path and the same answer
      self.assertEqual(s0.stats()["iter_count"], s1.stats()["iter_count"])
      for k, d in self.NXC_DIGITS:
        self.checkarray(r0[k], r1[k], sd+":nxc-vs-plain:"+k, digits=d)

      # against ipopt; ipmc relaxes bounds by 1e-8, ipopt does not: digits 4 to 6
      for k, d in [("f", 4), ("x", 5), ("g", 6), ("lam_g", 4), ("lam_x", 4)]:
        self.checkarray(r_ref[k], r1[k], sd+":nxc-vs-ipopt:"+k, digits=d)

  @requires_nlpsol("ipmc")
  def test_ipmc_nxc_horizons(self):
    self.message("ipmc nxc on short horizons and partial constant blocks")
    # N=1 has a single transition, both first and last step of the recursion
    for N, nc in [(1, 3), (2, 2), (10, 4)]:
      p = self.nxc_ocp(N=N, nc=nc, p_goal=min(2.0, 0.25*N))
      for sd in ["manual", "auto"]:
        s0, r0 = self.nxc_solve(p, sd, None)
        self.assertTrue(s0.stats()["success"])
        for nxc in range(nc+1):
          tag = "N=%d nc=%d %s nxc=%d" % (N, nc, sd, nxc)
          s1, r1 = self.nxc_solve(p, sd, nxc)
          self.assertTrue(s1.stats()["success"], tag)
          self.assertEqual(s0.stats()["iter_count"], s1.stats()["iter_count"], tag)
          if nxc == 0:
            for k in ["f", "x", "g", "lam_g", "lam_x"]:
              self.assertEqual(float(ca.norm_inf(r0[k]-r1[k])), 0.0,
                               tag+": nxc=0 is not bit-identical on "+k)
          # the multipliers of these short horizons agree to ~1e-10, not 1e-11
          for k, d in self.NXC_DIGITS:
            if k in ("lam_g", "lam_x"): d = 9
            self.checkarray(r0[k], r1[k], tag+":"+k, digits=d)

  @requires_nlpsol("ipmc")
  def test_ipmc_nxc_paths(self):
    self.message("ipmc nxc through both structured branches and the delta_c factorization")
    for N, nc in [(5, 3), (10, 4)]:
      # terminal_eq pulls equality rows back through the dynamics
      for te in [False, True]:
        p = self.nxc_ocp(N=N, nc=nc, terminal_eq=te)
        # linsol_perturbed_mode uses the delta_c factorization
        for pert in [False, True]:
          extra = {"linsol_perturbed_mode": True} if pert else {}
          s0, r0 = self.nxc_solve(p, "auto", None, {"ipmc": dict(self.IPMC_OPTS, **extra)})
          s1, r1 = self.nxc_solve(p, "auto", p["nxc"],
                                  {"ipmc": dict(self.IPMC_OPTS, **extra)})
          tag = "N=%d nc=%d terminal_eq=%s perturbed=%s" % (N, nc, te, pert)
          self.assertTrue(s0.stats()["success"], tag)
          self.assertTrue(s1.stats()["success"], tag)
          self.assertEqual(s0.stats()["iter_count"], s1.stats()["iter_count"], tag)
          for k, d in self.NXC_DIGITS:
            self.checkarray(r0[k], r1[k], tag+":"+k, digits=d)

  @requires_nlpsol("ipmc")
  def test_ipmc_nxc_rejected(self):
    self.message("ipmc refuses an nxc the problem does not satisfy")
    p = self.nxc_ocp(N=6, nc=3)
    N, nc = p["N"], p["nc"]

    # no horizon at structure_detection none, so even nxc=0 is refused
    with self.assertInException("structure_detection"):
      ca.nlpsol("solver", "ipmc", p["nlp"],
                {"structure_detection": "none", "nxc": 0, "ipmc": dict(self.IPMC_OPTS)})
    # nxc is one number for the whole horizon, not one per stage
    with self.assertInException("cannot be cast to OT_INT"):
      self.nxc_solve(p, "auto", [nc]*(N+1))
    # more constant states than states (nx is 2+nc)
    with self.assertInException("is not in [0, nx["):
      self.nxc_solve(p, "auto", 3+nc)
    # negative
    with self.assertInException("is not in [0, nx["):
      self.nxc_solve(p, "auto", -1)

    # nxc=nc+1 declares v constant, whose dynamics are not the identity
    with self.assertInException("declared constant"):
      self.nxc_solve(p, "auto", nc+1)
    # the generated code refuses it as well
    if IPMC_CODEGEN:
      solver = ca.nlpsol("ipmc_nxc_rejected", "ipmc", p["nlp"],
                         {"structure_detection": "auto", "equality": p["equality"],
                          "nxc": nc+1, "ipmc": dict(self.IPMC_OPTS)})
      solver.generate("ipmc_nxc_rejected.c")
      r = self.compile_external("ipmc_nxc_rejected", "ipmc_nxc_rejected.c", **IPMC_CODEGEN)
      if r is not None:
        with self.assertRaises(Exception):
          r[0](**p["bounds"])
        os.remove(r[1])
      os.remove("ipmc_nxc_rejected.c")

  @requires_nlpsol("ipmc")
  def test_ipmc_nxc_restoration(self):
    self.message("ipmc nxc through the restoration phase")
    # p_goal=20 with |u| <= 1 is unreachable, so ipmc ends in restoration
    for N, nc, te in [(10, 4, True), (10, 4, False)]:
      p = self.nxc_ocp(N=N, nc=nc, p_goal=20.0, u_max=1.0, terminal_eq=te)
      for sd in ["auto", "manual"]:
        s0, r0 = self.nxc_solve(p, sd, None)
        st0 = s0.stats()
        tag0 = "N=%d nc=%d terminal_eq=%s %s" % (N, nc, te, sd)
        # restoration is reached
        self.assertTrue(st0["ipmc"]["restoration_iterations_count"] > 0,
                        tag0+": no restoration iterations, fixture is not "
                             "testing what it says it tests")
        self.assertFalse(st0["success"], tag0)
        for nxc in range(nc+1):
          tag = tag0 + " nxc=%d" % nxc
          s1, r1 = self.nxc_solve(p, sd, nxc)
          st1 = s1.stats()
          # the same iterate path, restoration included
          self.assertEqual(st0["iter_count"], st1["iter_count"], tag)
          self.assertEqual(st0["ipmc"]["restoration_iterations_count"],
                           st1["ipmc"]["restoration_iterations_count"], tag)
          self.assertEqual(st0["ipmc"]["return_flag"],
                           st1["ipmc"]["return_flag"], tag)
          if nxc == 0:
            for k in ["f", "x", "g", "lam_g", "lam_x"]:
              self.assertEqual(float(ca.norm_inf(r0[k]-r1[k])), 0.0,
                               tag+": nxc=0 is not bit-identical on "+k)
          for k, d in self.NXC_RESTO_DIGITS:
            self.checkarray(r0[k], r1[k], tag+":"+k, digits=d)

  @requires_nlpsol("ipmc")
  def test_ipmc_nxc_codegen(self):
    self.message("ipmc nxc through codegen and serialization")
    p = self.nxc_ocp(N=5, nc=3, p_goal=1.5)
    solver, r = self.nxc_solve(p, "auto", p["nxc"])
    self.assertTrue(solver.stats()["success"])
    if IPMC_CODEGEN:
      self.check_codegen(solver, p["bounds"], **IPMC_CODEGEN)
    self.check_serialize(solver, p["bounds"])

  @requires_nlpsol("ipmc")
  def test_ipmc_nxc_resolve(self):
    self.message("ipmc nxc on one solver object with changing bounds")
    p = self.nxc_ocp(N=8, nc=3)
    solver, _ = self.nxc_solve(p, "auto", p["nxc"])
    plain = ca.nlpsol("plain", "ipmc", p["nlp"],
                      {"structure_detection": "auto", "print_time": False,
                       "equality": p["equality"], "ipmc": dict(self.IPMC_OPTS)})
    for goal in [1.0, 1.5, 2.0]:
      b = dict(p["bounds"])
      lbg = list(b["lbg"]); lbg[-1] = goal; b["lbg"] = lbg
      r1 = solver(**b)
      r0 = plain(**b)
      self.assertTrue(solver.stats()["success"])
      self.assertEqual(solver.stats()["iter_count"], plain.stats()["iter_count"])
      for k, d in self.NXC_DIGITS:
        self.checkarray(r0[k], r1[k], "goal=%g:%s" % (goal, k), digits=d)

  # soft constraints (nlpsol S, s and f_s) on a structured OCP
  IPMC_SLACK_OPTS = {"print_level": 0, "tol": 1e-9, "max_iter": 500}

  # On a capped or L-infinity slack the optimum is nearly flat
  # in some directions and ipopt may stop at an "acceptable" point ~1e-3 off
  # in x.  The reference takes the first of these variants that reports
  # Solve_Succeeded: the adaptive barrier update, or a start close to the
  # bounds, gets through where the default stalls.
  IPOPT_REF_VARIANTS = [{}, {"mu_strategy": "adaptive"},
                        {"bound_push": 1e-6, "bound_frac": 1e-6}]
  # ipopt reference on the problem as written, without bound relaxation
  IPOPT_REF_OPTS = {"print_time": False,
                    "ipopt": {"print_level": 0, "sb": "yes", "tol": 1e-12,
                              "constr_viol_tol": 1e-12, "dual_inf_tol": 1e-9,
                              "compl_inf_tol": 1e-12, "acceptable_tol": 1e-12,
                              "bound_relax_factor": 0, "max_iter": 3000}}

  def slack_ocp(self, N, soften="band", group_mode="per_row", layout="pair",
                v_lo=1.0, v_hi=7.0, p_goal=5.0, T=1.0, u_max=40.0,
                slope=5.5) -> Dict[str, Any]:
    """Sprint OCP with softened rows; soften picks the rows, group_mode and layout S.

    S = [S_lo; S_up] is built by slack_layout from the softened rows, each
    tagged with its stage and its kind ("band": the speed band, "cap": the
    position corridor)."""
    dt = T/N
    xs = ca.SX.sym("xs", 2)
    us = ca.SX.sym("us")
    ode = ca.Function("ode", [xs, us], [ca.vertcat(xs[1], us-0.02*xs[1]*xs[1])])
    k1 = ode(xs, us)
    k2 = ode(xs+dt/2*k1, us)
    k3 = ode(xs+dt/2*k2, us)
    k4 = ode(xs+dt*k3, us)
    F = ca.Function("F", [xs, us], [xs+dt/6*(k1+2*k2+2*k3+k4)])

    X = [ca.SX.sym("x_%d" % k, 2) for k in range(N+1)]
    U = [ca.SX.sym("u_%d" % k) for k in range(N)]

    lbx, ubx = [], []
    for k in range(N+1):
      lbx += [0.0, 0.0] if k == 0 else [-inf, -inf]
      ubx += [0.0, 0.0] if k == 0 else [inf, inf]
      if k < N:
        lbx.append(-u_max)
        ubx.append(u_max)

    var, g, lbg, ubg, ng, equality = [], [], [], [], [], []
    soft = []   # (row index in the stacked [g; x], stage, kind)

    # corridor: two path rows per stage sharing one column
    def paths(k):
      rows = []
      if soften == "corridor":
        rows = [(X[k][1], -inf, v_hi, "band"), (X[k][0], -inf, slope*k*dt, "cap")]
      elif soften == "band_cap":
        rows = [(X[k][1], v_lo if k >= 1 else -inf, v_hi, "band"),
                (X[k][0], -inf, slope*k*dt, "cap")]
      elif soften == "band":
        rows = [(X[k][1], v_lo if k >= 1 else -inf, v_hi, "band")]
      elif soften == "upper":
        rows = [(X[k][1], -inf, v_hi, "band")]
      n = 0
      for expr, lb, ub, kind in rows:
        if k >= 1:
          soft.append((len(equality), k, kind))
        g.append(expr)
        lbg.append(lb)
        ubg.append(ub)
        equality.append(lb == ub)
        n += 1
      if k == N:   # terminal reach constraint, always hard
        g.append(X[N][0])
        lbg.append(p_goal)
        ubg.append(inf)
        equality.append(False)
        n += 1
      return n

    for k in range(N):
      var += [X[k], U[k]]
      g.append(X[k+1]-F(X[k], U[k]))     # gap-closing / dynamics rows
      lbg += [0.0, 0.0]
      ubg += [0.0, 0.0]
      equality += [True, True]
      ng.append(paths(k))
    var.append(X[N])
    ng.append(paths(N))

    x = ca.vertcat(*var)
    G = ca.vertcat(*g)
    nxu, ngu = x.numel(), G.numel()

    # bound_x and bound_band_x soften a one- or two-sided simple bound on v
    if soften in ["bound_x", "bound_band_x"]:
      for k in range(1, N+1):
        ubx[3*k+1] = v_hi
        if soften == "bound_band_x":
          lbx[3*k+1] = v_lo
        soft.append((ngu+3*k+1, k, "band"))

    S, lo_cols, up_cols = self.slack_layout(ngu+nxu, soft, layout, group_mode)
    n_lift = self.slack_n_lift

    return dict(x=x, f=dt*sum(U[k]*U[k] for k in range(N)), g=G,
                nxu=nxu, ngu=ngu, ns=S.size2(), S=S, lo_cols=lo_cols, up_cols=up_cols,
                n_lift=n_lift,
                lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg, equality=equality,
                structure=dict(N=N, nx=[2]*(N+1), nu=[1]*N+[0], ng=ng))

  # How the columns of S lay over the softened rows.  A layout maps (kind,
  # side, stage) to a column key, or None for a hard side; equal keys are one
  # column.  g is the group_mode key of a stage: "single" one budget for all
  # stages (L-infinity, lifted into a helper state once the stages are detected),
  # "mixed" one budget for the odd stages and one column per even stage,
  # "per_row" one column per stage.
  #   pair        a column per side: [S 0; 0 S] (s = [s_l; s_u])
  #   sym         one column for both sides of a row: [S; S]
  #   lo, up      one side only: [S; 0], [0; S]
  #   stage_mixed one column per stage on the LOWER side of its speed band and
  #               the UPPER side of its position cap: shared within a stage,
  #               on mixed sides
  #   linf_mixed  the same, one column for all stages: lifted, mixed sides
  #   glob_loc    the lower side of every band row by a column of its stage,
  #               the upper side by one lifted column: a helper and a
  #               stage-local slack on the same row
  slack_layouts = {
    "pair":        lambda kind, side, g: (side, g),
    "sym":         lambda kind, side, g: ("s", g),
    "lo":          lambda kind, side, g: ("lo", g) if side == "lo" else None,
    "up":          lambda kind, side, g: ("up", g) if side == "up" else None,
    "stage_mixed": lambda kind, side, g: ("m", g) if (kind, side) in
                   [("band", "lo"), ("cap", "up")] else None,
    "linf_mixed":  lambda kind, side, g: ("m", "all") if (kind, side) in
                   [("band", "lo"), ("cap", "up")] else None,
    "glob_loc":    lambda kind, side, g: (("lo", g) if side == "lo" else ("up", "all"))
                   if kind == "band" else None,
  }

  def slack_layout(self, n, soft, layout, group_mode="per_row"):
    """S = [S_lo; S_up] (2n rows), and the columns on lower and upper sides"""
    def group(k):
      if group_mode == "single":
        return "all"
      if group_mode == "mixed":
        return "all" if k % 2 else ("row", k)
      return ("row", k)
    f = self.slack_layouts[layout]
    columns, rows_, cols_, lo_cols, up_cols = {}, [], [], set(), set()
    stages = {}
    # all lower sides first, so that a pair layout numbers [s_l; s_u]
    for side in ["lo", "up"]:
      for r, k, kind in soft:
        key = f(kind, side, group(k))
        if key is None: continue
        if key not in columns:
          columns[key] = len(columns)
        rows_.append(r if side == "lo" else n+r)
        cols_.append(columns[key])
        (lo_cols if side == "lo" else up_cols).add(columns[key])
        stages.setdefault(columns[key], set()).add(k)
    S = ca.Sparsity.triplet(2*n, len(columns), rows_, cols_)
    # columns whose rows span several stages: lifted once stages are detected
    self.slack_n_lift = sum(1 for v in stages.values() if len(v) > 1)
    return S, sorted(lo_cols), sorted(up_cols)

  def slack_case(self, penalty, w, ubs=None, par=None,
                 **kwargs) -> Dict[str, Any]:
    """slack_ocp() plus the penalty f_s, which may depend on the parameter par"""
    p = self.slack_ocp(**kwargs)
    s = ca.SX.sym("s", p["ns"])
    p["s"] = s
    p["ubs"] = ubs
    p["nlp"] = {"x": p["x"], "f": p["f"], "g": p["g"], "s": s, "f_s": penalty(s, w)}
    p["hard_nlp"] = {"x": p["x"], "f": p["f"], "g": p["g"]}
    p["bounds"] = dict(x0=0, lbx=p["lbx"], ubx=p["ubx"],
                       lbg=p["lbg"], ubg=p["ubg"])
    if par is not None:
      p["par"] = par
      p["nlp"]["p"] = par
    return p

  def slack_activity(self, s, p, tol=1e-3):
    """Number of active columns on lower and on upper sides; tol is far above
    the L2 barrier floor 5e-6"""
    return (sum(1 for j in p["lo_cols"] if float(s[j]) > tol),
            sum(1 for j in p["up_cols"] if float(s[j]) > tol))

  def ipmc_slack_solver(self, p, structure_detection):
    opts = {"structure_detection": structure_detection,
            "ipmc": dict(self.IPMC_SLACK_OPTS),
            "S": p["S"]}
    if structure_detection == "manual":
      opts.update(p["structure"])
    elif structure_detection == "auto":
      opts["equality"] = p["equality"]
    return ca.nlpsol("solver", "ipmc", p["nlp"], opts)

  def ipmc_hard_solve(self, p):
    """The hard twin, dense. Returns (result, success)."""
    solver = ca.nlpsol("hard", "ipmc", p["hard_nlp"],
                       {"structure_detection": "none", "ipmc": dict(self.IPMC_SLACK_OPTS)})
    r = solver(**p["bounds"])
    return r, solver.stats()["success"]

  def slack_ocp_reference(self, p, pval=None):
    """Augmentation relaxing every row, solved dense by ipopt"""
    x, g, s = p["x"], p["g"], p["s"]
    nx, ng, ns = p["nxu"], p["ngu"], p["ns"]
    Sd = ca.DM.ones(p["S"])
    n = nx+ng
    S_lo, S_up = Sd[:n, :], Sd[n:, :]
    G = ca.vertcat(g+ca.mtimes(S_lo[:ng, :], s),    # >= lbg
                   g-ca.mtimes(S_up[:ng, :], s),    # <= ubg
                   x+ca.mtimes(S_lo[ng:, :], s),    # >= lbx
                   x-ca.mtimes(S_up[ng:, :], s))    # <= ubx
    lbG = ca.vertcat(ca.DM(p["lbg"]), -inf*ca.DM.ones(ng),
                     ca.DM(p["lbx"]), -inf*ca.DM.ones(nx))
    ubG = ca.vertcat(inf*ca.DM.ones(ng), ca.DM(p["ubg"]),
                     inf*ca.DM.ones(nx), ca.DM(p["ubx"]))
    ubs = inf*ca.DM.ones(ns) if p["ubs"] is None else ca.DM(p["ubs"])

    rnlp = {"x": ca.vertcat(x, s), "f": p["f"]+p["nlp"]["f_s"], "g": G}
    extra = {}
    if "par" in p:
      rnlp["p"] = p["par"]
      extra["p"] = pval
    args = dict(lbx=ca.vertcat(-inf*ca.DM.ones(nx), ca.DM.zeros(ns)),
                ubx=ca.vertcat(inf*ca.DM.ones(nx), ubs), lbg=lbG, ubg=ubG, **extra)
    for variant in self.IPOPT_REF_VARIANTS:
      solver = ca.nlpsol("reference", "ipopt", rnlp,
                         dict(self.IPOPT_REF_OPTS,
                              ipopt=dict(self.IPOPT_REF_OPTS["ipopt"], **variant)))
      r = solver(x0=0, **args)
      if solver.stats()["return_status"] == "Solve_Succeeded": break
    self.assertEqual(solver.stats()["return_status"], "Solve_Succeeded")
    return {"f": r["f"], "x": r["x"][:nx], "s": r["x"][nx:],
            "g": ca.Function("g", [x], [g])(r["x"][:nx]),
            "lam_s": r["lam_x"][nx:]}

  def slack_feasible(self, p, r, tol=1e-6):
    """r is feasible for the relaxed bounds, hard rows included"""
    Sd = ca.DM.ones(p["S"])
    n = p["nxu"]+p["ngu"]
    z = ca.vertcat(r["g"], r["x"])
    lb = ca.vertcat(ca.DM(p["lbg"]), ca.DM(p["lbx"]))
    ub = ca.vertcat(ca.DM(p["ubg"]), ca.DM(p["ubx"]))
    self.assertTrue(float(ca.mmax(z-ub-ca.mtimes(Sd[n:, :], r["s"]))) < tol)
    self.assertTrue(float(ca.mmax(lb-ca.mtimes(Sd[:n, :], r["s"])-z)) < tol)

  def check_slack_structure(self, solver, p, structure_detection):
    """The detected OCP structure is the user's, slacks aside"""
    if structure_detection == "none":
      return
    st = solver.stats()
    self.assertEqual(st["N"], p["structure"]["N"])
    for k in ["nx", "nu", "ng"]:
      self.checkarray(ca.DM(p["structure"][k]), ca.DM(st[k]),
                      "structure:"+k, digits=12)

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  @memory_heavy()
  def test_ipmc_slacks(self):
    self.message("ipmc slacks: L1, L2, softened simple bounds, a within-stage shared column")
    L1 = lambda s, w: w*ca.sum1(s)
    L2 = lambda s, w: w*ca.dot(s, s)
    cases = [
      # name, penalty, w, kwargs, (ns, nnz(S)), min df, min max|s|, min active (lo, up)
      ("l1_band_N20",     L1, 0.5, dict(N=20), (40, 40), 1.0, 0.1, (1, 5)),
      ("l2_band_N20",     L2, 2.0, dict(N=20), (40, 40), 1.0, 0.1, (1, 5)),
      ("l1_bound_x_N10",  L1, 0.5, dict(N=10, soften="bound_x", v_hi=6.0),
       (20, 20), 10.0, 1.0, (0, 5)),
      # a two-sided simple bound: both sides are active at the optimum
      ("l1_bound_band_x_N10", L1, 0.5,
       dict(N=10, soften="bound_band_x", v_lo=3.0, v_hi=6.0),
       (20, 20), 10.0, 1.0, (2, 5)),
      # the same, one symmetric column per bound
      ("l1_bound_band_x_sym_N10", L1, 0.5,
       dict(N=10, soften="bound_band_x", v_lo=3.0, v_hi=6.0, layout="sym"),
       (10, 20), 10.0, 1.0, (7, 7)),
      # two path rows per stage sharing one column, which never leaves the stage
      ("stage_shared_N8", L1, 1.0, dict(N=8, soften="corridor", v_hi=6.0),
       (16, 32), 10.0, 1.0, (0, 4)),
    ]
    for name, penalty, w, kwargs, shape, df_min, s_min, act in cases:
      p = self.slack_case(penalty, w, **kwargs)
      self.assertEqual((p["ns"], p["S"].nnz()), shape)
      ref = self.slack_ocp_reference(p)
      rh, hard_ok = self.ipmc_hard_solve(p)
      self.assertTrue(hard_ok)

      # softening pays off, and both sides are genuinely active
      nsl, nsu = self.slack_activity(ref["s"], p)
      print("test_ipmc_slacks", name, "f_soft", float(ref["f"]),
            "f_hard", float(rh["f"]), "max|s|", float(ca.norm_inf(ref["s"])),
            "active", (nsl, nsu))
      self.assertTrue(float(rh["f"]-ref["f"]) > df_min)
      self.assertTrue(float(ca.norm_inf(ref["s"])) > s_min)
      self.assertTrue(nsl >= act[0] and nsu >= act[1])

      for sd in ["none", "manual", "auto"]:
        print("test_ipmc_slacks", name, sd)
        solver = self.ipmc_slack_solver(p, sd)
        r = solver(**p["bounds"])
        self.assertTrue(solver.stats()["success"])
        self.check_slack_structure(solver, p, sd)
        # ipmc relaxes bounds by 1e-8, ipopt does not; L2 slacks stall at 8.8e-6: s digits 4
        for k, d in [("f", 4), ("x", 5), ("g", 6), ("s", 4)]:
          self.checkarray(ref[k], r[k], name+":"+sd+":"+k, digits=d)
        # every column is local to its stage
        self.assertEqual(solver.stats()["n_lift"], 0)
        # feasible for the relaxed bounds, hard rows included
        self.slack_feasible(p, r)
        # z and Z are values in the serialized stream
        self.check_serialize(solver, p["bounds"])

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  @memory_heavy()
  def test_ipmc_slacks_layouts(self):
    self.message("ipmc slacks: every column layout x penalty x cap x structure_detection")
    # (name, slack_ocp kwargs); see slack_layouts
    layouts = [
      ("pair",          dict(layout="pair")),
      ("sym",           dict(layout="sym")),
      ("lo",            dict(layout="lo")),
      ("up",            dict(layout="up")),
      ("stage_mixed",   dict(layout="stage_mixed", soften="band_cap", slope=4.8)),
      # lifted: one budget for both sides, two budgets, mixed sides, a helper
      # and a local slack on one row, helpers next to locals
      ("linf_sym",      dict(layout="sym", group_mode="single")),
      ("linf_two",      dict(layout="pair", group_mode="single")),
      ("linf_up",       dict(layout="up", group_mode="single")),
      ("linf_mixed",    dict(layout="linf_mixed", soften="band_cap", slope=4.8)),
      ("glob_loc",      dict(layout="glob_loc")),
      ("linf_and_local", dict(layout="sym", group_mode="mixed")),
    ]
    penalties = [
      ("L1",   lambda s, w: 0.5*ca.sum1(s)),
      ("L2",   lambda s, w: 2.0*ca.dot(s, s)),
      ("L1L2", lambda s, w: ca.dot(ca.DM([0.3+0.1*j for j in range(s.numel())]), s)
                            + 0.5*ca.dot(s, s)),
    ]
    for lname, lkw in layouts:
      checked = False
      for pname, penalty in penalties:
        for cap in [None, "finite"]:
          if cap is not None and pname != "L1": continue
          kw = dict(dict(N=10, v_lo=3.0, v_hi=6.0), **lkw)
          p = self.slack_case(penalty, None, **kw)
          ns = p["ns"]
          ubs = None if cap is None else ca.DM([0.2+0.05*j for j in range(ns)])
          p["ubs"] = ubs
          tag = "%s/%s%s" % (lname, pname, "" if cap is None else "/ubs")
          ref = self.slack_ocp_reference(p)
          nlo, nup = self.slack_activity(ref["s"], p)
          print("test_ipmc_slacks_layouts", tag, "ns", ns, "lifted", p["n_lift"],
                "active", (nlo, nup), "f", float(ref["f"]))
          # something is relaxed, or the case proves nothing
          self.assertTrue(nlo+nup >= 1, tag)
          if lname.startswith("linf") or lname == "glob_loc":
            self.assertTrue(p["n_lift"] >= 1, tag)
          args = dict(p["bounds"])
          if ubs is not None: args["ubs"] = ubs
          for sd in ["none", "manual", "auto"]:
            solver = self.ipmc_slack_solver(p, sd)
            r = solver(**args)
            self.assertTrue(solver.stats()["success"], tag+":"+sd)
            self.check_slack_structure(solver, p, sd)
            self.assertEqual(solver.stats()["n_lift"],
                             0 if sd == "none" else p["n_lift"], tag+":"+sd)
            for k, d in [("f", 4), ("x", 5), ("g", 5), ("s", 4)]:
              self.checkarray(ref[k], r[k], tag+":"+sd+":"+k, digits=d)
            self.slack_feasible(p, r)
            if ubs is not None:
              self.assertTrue(float(ca.mmax(r["s"]-ubs)) < 1e-6, tag)
            if sd == "auto" and cap is None and pname == "L1":
              # a starting slack changes the route, not the answer
              r0 = solver(**dict(args, s0=0.5*ca.DM.ones(ns)))
              self.assertTrue(solver.stats()["success"], tag)
              for k, d in [("f", 5), ("x", 5), ("s", 4)]:
                self.checkarray(r[k], r0[k], tag+":s0:"+k, digits=d)
            # codegen and serialization once per layout
            if sd == "auto" and not checked:
              checked = True
              if IPMC_CODEGEN:
                self.check_codegen(solver, args, **IPMC_CODEGEN)
              self.check_serialize(solver, args)

  @requires_nlpsol("ipmc")
  def test_ipmc_slacks_degenerate(self):
    self.message("ipmc slacks: ubs = 0 is refused, inactive slacks give the hard answer")
    # ground truth is the hard NLP, so this test runs without ipopt
    L1 = lambda s, w: w*ca.sum1(s)
    p = self.slack_case(L1, 0.5, ubs=ca.DM.zeros(2*20), N=20)
    self.assertEqual(p["ns"], 2*20)
    for sd in ["none", "manual", "auto"]:
      solver = self.ipmc_slack_solver(p, sd)
      with self.assertInException("ubs = 0"):
        solver(**dict(p["bounds"], ubs=p["ubs"]))
    cases = [
      ("inactive_N10", dict(N=10, v_lo=-50.0, v_hi=50.0), None),
    ]
    for name, kwargs, ubs in cases:
      p = self.slack_case(L1, 0.5, ubs=ubs, **kwargs)
      self.assertEqual(p["ns"], 2*kwargs["N"])
      rh, hard_ok = self.ipmc_hard_solve(p)
      self.assertTrue(hard_ok)
      for sd in ["none", "manual", "auto"]:
        print("test_ipmc_slacks_degenerate", name, sd)
        solver = self.ipmc_slack_solver(p, sd)
        args = dict(p["bounds"])
        if ubs is not None:
          args["ubs"] = ubs
        r = solver(**args)
        self.assertTrue(solver.stats()["success"])
        self.check_slack_structure(solver, p, sd)
        # the slacks are inert
        self.checkarray(ca.DM.zeros(p["ns"]), r["s"], name+":"+sd+":s", digits=6)
        # so the answer is the hard one, up to the barrier floor (1.8e-9)
        for k in ["f", "x", "g"]:
          self.checkarray(rh[k], r[k], name+":"+sd+":"+k, digits=7)

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  def test_ipmc_slacks_infeasible(self):
    self.message("ipmc slacks make an infeasible OCP solvable")
    # p(T)=5 needs a mean speed of 5 but the cap is 4, so the hard twin fails
    p = self.slack_case(lambda s, w: w*ca.sum1(s), 1.0,
                        N=25, soften="upper", v_hi=4.0, p_goal=5.0)
    rh, hard_ok = self.ipmc_hard_solve(p)
    self.assertFalse(hard_ok)
    ref = self.slack_ocp_reference(p)
    self.assertTrue(float(ca.norm_inf(ref["s"])) > 1.0)

    for sd in ["none", "manual", "auto"]:
      print("test_ipmc_slacks_infeasible", sd)
      solver = self.ipmc_slack_solver(p, sd)
      r = solver(**p["bounds"])
      self.assertTrue(solver.stats()["success"])
      self.check_slack_structure(solver, p, sd)
      for k, d in [("f", 4), ("x", 5), ("g", 6), ("s", 5)]:
        self.checkarray(ref[k], r[k], "infeasible_hard_N25:"+sd+":"+k, digits=d)

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  @memory_heavy()
  def test_ipmc_slacks_shared_column(self):
    self.message("ipmc slacks: one Linf budget per side shared by every stage")
    # at manual and auto the columns are lifted into constant helper states
    L1 = lambda s, w: w*ca.sum1(s)
    cases = [
      ("linf_band_N20", 1.0, dict(N=20, group_mode="single"), 1.0, 0.1, (1, 1)),
      ("linf_infeasible_N12", 5.0,
       dict(N=12, group_mode="single", soften="upper", v_hi=4.0, p_goal=5.0),
       None, 1.0, (0, 1)),
    ]
    for name, w, kwargs, df_min, s_min, act in cases:
      p = self.slack_case(L1, w, **kwargs)
      self.assertEqual(p["ns"], 2)
      ref = self.slack_ocp_reference(p)
      rh, hard_ok = self.ipmc_hard_solve(p)
      nsl, nsu = self.slack_activity(ref["s"], p)
      self.assertTrue(float(ca.norm_inf(ref["s"])) > s_min)
      self.assertTrue(nsl >= act[0] and nsu >= act[1])
      if df_min is None:
        self.assertFalse(hard_ok)   # infeasible hard twin
      else:
        self.assertTrue(hard_ok)
        self.assertTrue(float(rh["f"]-ref["f"]) > df_min)

      for sd in ["none", "manual", "auto"]:
        print("test_ipmc_slacks_shared_column", name, sd, (nsl, nsu))
        solver = self.ipmc_slack_solver(p, sd)
        r = solver(**p["bounds"])
        self.assertTrue(solver.stats()["success"])
        self.check_slack_structure(solver, p, sd)
        # nothing is lifted at none
        self.assertEqual(solver.stats()["n_lift"], 0 if sd == "none" else 2)
        for k, d in [("f", 4), ("x", 5), ("g", 6), ("s", 5)]:
          self.checkarray(ref[k], r[k], name+":"+sd+":"+k, digits=d)
        # the lifted description survives serialization
        self.check_serialize(solver, p["bounds"])

  @requires_nlpsol("ipmc")
  def test_ipmc_unheld_multipliers_zero(self):
    self.message("ipmc: unbounded entries have multiplier 0, independent of the dual guess")
    x = ca.SX.sym("x", 3)
    nlp = {"x": x, "f": (x[0]-1)**2 + (x[1]-2)**2 + x[2]**2 + x[0]*x[2],
           "g": ca.vertcat(x[0]+x[1], x[0]*x[2])}
    bounds = dict(x0=0.5, lbx=[-5, -inf, -inf], ubx=[5, inf, inf],
                  lbg=[-10, -inf], ubg=[10, inf])
    guess = dict(lam_x0=[0.3, -0.7, 0.9], lam_g0=[0.4, -1.3])
    solver = ca.nlpsol("solver", "ipmc", nlp, {"ipmc": {"print_level": 0}})
    cases = [("plain", solver, bounds, [1, 2], [1])]

    L1 = lambda s, w: w*ca.sum1(s)
    q = self.slack_case(L1, 1.0, N=20, group_mode="single")
    lifted = self.ipmc_slack_solver(q, "auto")
    lbx, ubx = numpy.array(ca.DM(q["lbx"])).ravel(), numpy.array(ca.DM(q["ubx"])).ravel()
    free_x = [i for i in range(len(lbx)) if lbx[i] == -inf and ubx[i] == inf]
    self.assertTrue(len(free_x) > 0)
    cases.append(("lifted", lifted, q["bounds"], free_x, []))

    for tag, sv, b, free_x, free_g in cases:
      rng = numpy.random.default_rng(3)
      g = guess if tag == "plain" else dict(
        lam_x0=rng.standard_normal(sv.nnz_in("lam_x0")),
        lam_g0=rng.standard_normal(sv.nnz_in("lam_g0")))
      r0 = sv(**b)
      r = sv(**dict(b, **g))
      self.assertTrue(sv.stats()["success"], tag)
      if tag == "lifted":
        self.assertEqual(sv.stats()["n_lift"], 2)
      for i in free_x:
        self.assertEqual(float(r["lam_x"][i]), 0.0, tag)
      for i in free_g:
        self.assertEqual(float(r["lam_g"][i]), 0.0, tag)
      for k in r0.keys():
        self.assertTrue(numpy.array_equal(numpy.array(r0[k]), numpy.array(r[k])), tag+":"+k)
      if IPMC_CODEGEN:
        self.check_codegen(sv, dict(b, **g), **IPMC_CODEGEN)
      self.check_serialize(sv, dict(b, **g))

  def ipmc_small_ocp(self, N=3, free_row=False):
    """x_{k+1} = x_k + u_k, |u_k| <= 1, x_k + u_k <= 2, x_0 = 0; free_row adds free rows"""
    X = [ca.SX.sym("x%d" % k) for k in range(N+1)]
    U = [ca.SX.sym("u%d" % k) for k in range(N)]
    w, lbw, ubw, g, lbg, ubg, eq = [], [], [], [], [], [], []
    for k in range(N+1):
      w += [X[k]]; lbw += [0.0 if k == 0 else -inf]; ubw += [0.0 if k == 0 else inf]
      if k < N:
        w += [U[k]]; lbw += [-1.0]; ubw += [1.0]
        g += [X[k+1] - (X[k] + U[k])]; lbg += [0.0]; ubg += [0.0]; eq += [True]
        g += [X[k] + U[k]]; lbg += [-inf]; ubg += [2.0]; eq += [False]
        if free_row:
          g += [ca.sin(X[k])*U[k] + X[k]**3]; lbg += [-inf]; ubg += [inf]; eq += [False]
    f = sum((X[k]-1)**2 for k in range(N+1)) + sum(u**2 for u in U)
    return ({"x": ca.vertcat(*w), "f": f, "g": ca.vertcat(*g)},
            dict(lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg), eq)

  @requires_nlpsol("ipmc")
  def test_ipmc_equality_checks(self):
    self.message("ipmc: the equality option finds the gap-closing rows, the bounds the rest")
    nlp, b, eq = self.ipmc_small_ocp()
    def mk(e, extra={}):
      return ca.nlpsol("s", "ipmc", nlp, dict({"structure_detection": "auto", "equality": e,
                                               "print_time": False, "ipmc": dict(self.IPMC_OPTS)},
                                              **extra))
    with self.assertInException("requires the 'equality' option"):
      ca.nlpsol("s", "ipmc", nlp, {"structure_detection": "auto"})
    with self.assertInException("Expected 6 elements"):
      mk(eq+[True])
    # the fixed initial state needs no tag; the solve succeeds, twice on one object
    solver = mk(eq)
    r = solver(**b)
    self.assertTrue(solver.stats()["success"])
    r2 = solver(**b)
    self.assertTrue(numpy.array_equal(numpy.array(r["x"]), numpy.array(r2["x"])))
    # an untagged path row with lbg == ubg is an equality for that solve
    r2 = solver(**dict(b, lbg=[0, 0.5, 0, -inf, 0, -inf], ubg=[0, 0.5, 0, 2, 0, 2]))
    self.assertTrue(solver.stats()["success"])
    self.checkarray(r2["g"][1], 0.5, "g[1]", digits=10)
    # a tag on a path row without equal bounds changes nothing
    e = list(eq); e[1] = True
    r2 = mk(e)(**b)
    for k in ["x", "f", "lam_x", "lam_g"]:
      self.assertTrue(numpy.array_equal(numpy.array(r[k]), numpy.array(r2[k])), k)
    # a gap-closing row
    with self.assertInException("constraint row g[0] closes a gap of the dynamics"):
      mk(eq)(**dict(b, ubg=[0.5]+b["ubg"][1:]))
    # slacks: ubs must be positive; one column per side of the path rows
    S = ca.Sparsity.triplet(2*(6+7), 6, [1, 3, 5, 13+1, 13+3, 13+5], range(6))
    s = ca.SX.sym("s", 6)
    soft = ca.nlpsol(
      "s", "ipmc", dict(nlp, s=s, f_s=ca.sum1(s)),
      {"structure_detection": "auto", "equality": eq, "S": S, "print_time": False,
       "ipmc": dict(self.IPMC_OPTS)})
    soft(**b)
    self.assertTrue(soft.stats()["success"])
    with self.assertInException("slack column 4 has the upper bound ubs = 0, but relaxes a finite bound"):
      soft(**dict(b, ubs=[inf, inf, inf, inf, 0, inf]))
    # generated C: the solve, and the same refusals, which return nonzero and name the entry
    if IPMC_CODEGEN:
      self.check_codegen(solver, b, **IPMC_CODEGEN)
      refusals = [
        (solver, dict(b, ubg=[0.5]+b["ubg"][1:]), "constraint row g[0] closes a gap of the dynamics"),
        (soft, dict(b, ubs=[inf, inf, inf, inf, 0, inf]),
         "slack column 4 has an upper bound ubs <= 0, but relaxes a finite bound")]
      for f, inputs, msg in refusals:
        self.check_codegen(f, inputs, main=True, main_return_code=[1], **IPMC_CODEGEN)
        if not args.run_slow: continue
        with open(f.name()+"_out.txt") as out:
          self.assertIn(msg, out.read())
    self.check_serialize(solver, b)

  @requires_nlpsol("ipmc")
  def test_ipmc_pinned_bounds(self):
    self.message("ipmc pins variables with lbx == ubx at solve time")
    nlp, b, eq = self.ipmc_small_ocp(N=5)
    nx = nlp["x"].numel()
    def mk():
      return ca.nlpsol("s", "ipmc", nlp, {"structure_detection": "auto", "equality": eq,
                                          "print_time": False, "ipmc": dict(self.IPMC_OPTS)})
    def pinned(pins):
      lbx, ubx = [-inf]*nx, [inf]*nx
      for i in range(1, nx, 2):
        lbx[i], ubx[i] = -1.0, 1.0
      for i, v in pins.items():
        lbx[i] = ubx[i] = v
      return dict(b, lbx=lbx, ubx=ubx)
    calls = [pinned({0: 0.0}), pinned({}), pinned({0: 0.5, nx-1: 1.0}), pinned({0: 0.0}),
             pinned({0: 0.0, 5: 0.3})]
    lag = ca.Function("lag", [nlp["x"]], [ca.gradient(nlp["f"], nlp["x"]),
                                          ca.jacobian(nlp["g"], nlp["x"]).T])
    # one object through varying pins gives what a fresh solver gives, bit for bit
    solver = mk()
    for i, args in enumerate(calls):
      r = solver(**args)
      rf = mk()(**args)
      for k in ["x", "f", "lam_x", "lam_g"]:
        self.assertTrue(numpy.array_equal(numpy.array(r[k]), numpy.array(rf[k])), "%d:%s" % (i, k))
      self.assertTrue(solver.stats()["success"], i)
      # the pinned multipliers close stationarity
      gf, jgt = lag(r["x"])
      stat = gf + ca.mtimes(jgt, r["lam_g"]) + r["lam_x"]
      self.checkarray(stat, ca.DM.zeros(nx), "%d:stationarity" % i, digits=8)
      for j in range(nx):
        if args["lbx"][j] == args["ubx"][j]:
          self.checkarray(r["x"][j], args["lbx"][j], "%d:x[%d]" % (i, j), digits=10)
    self.assertNotEqual(float(r["lam_x"][5]), 0.0)
    if IPMC_CODEGEN:
      for args in calls[2:]:
        self.check_codegen(solver, args, **IPMC_CODEGEN)
    self.check_serialize(solver, calls[2])

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  def test_ipmc_pinned_rows(self):
    self.message("ipmc holds path rows with lbg == ubg as equalities at solve time")
    nlp, b, eq = self.ipmc_small_ocp(N=5)
    def pinned(pins):
      lbg, ubg = list(b["lbg"]), list(b["ubg"])
      for k, v in pins.items():
        lbg[2*k+1] = ubg[2*k+1] = v
      return dict(b, lbg=lbg, ubg=ubg)
    calls = [pinned({2: 1.5}), pinned({}), pinned({1: 0.8, 4: 1.2}), pinned({2: 1.5}),
             pinned({2: 1.0})]
    ref = ca.nlpsol("reference", "ipopt", nlp,
                    dict(self.IPOPT_HARD_REF_OPTS, ipopt=dict(self.IPOPT_HARD_REF_OPTS["ipopt"])))
    for sd in ["none", "manual", "auto"]:
      opts = {"structure_detection": sd, "print_time": False, "ipmc": dict(self.IPMC_OPTS)}
      if sd == "manual":
        opts.update(N=5, nx=[1]*6, nu=[1]*5+[0], ng=[1]*5+[0])
      elif sd == "auto":
        opts["equality"] = eq
      # one object through varying pins gives what a fresh solver gives, bit for bit
      solver = ca.nlpsol("s", "ipmc", nlp, opts)
      for i, args in enumerate(calls):
        tag = "%s:%d" % (sd, i)
        r = solver(**args)
        self.assertTrue(solver.stats()["success"], tag)
        rf = ca.nlpsol("s", "ipmc", nlp, opts)(**args)
        for k in ["x", "f", "g", "lam_x", "lam_g"]:
          self.assertTrue(numpy.array_equal(numpy.array(r[k]), numpy.array(rf[k])), tag+":"+k)
        # ipmc relaxes bounds by 1e-8, ipopt does not
        r_ref = ref(**args)
        for k, d in [("f", 5), ("x", 5), ("g", 6), ("lam_g", 4), ("lam_x", 4)]:
          self.checkarray(r_ref[k], r[k], tag+":"+k, digits=d)
        # a held row pulls either way, which the row as an upper bound could not
        if i == 2:
          self.assertTrue(float(r["lam_g"][3]) > 0.1 and float(r["lam_g"][9]) < -0.1, tag)
    if IPMC_CODEGEN:
      for args in calls[2:]:
        self.check_codegen(solver, args, **IPMC_CODEGEN)
    self.check_serialize(solver, calls[2])

  @requires_nlpsol("ipmc")
  def test_ipmc_free_rows_inert(self):
    self.message("ipmc holds rows without bounds inert")
    # the free rows leave the iterate path unchanged, also truncated and perturbed
    nlp0, b0, eq0 = self.ipmc_small_ocp(N=5)
    nlp1, b1, eq1 = self.ipmc_small_ocp(N=5, free_row=True)
    free_rows = [3*k+2 for k in range(5)]
    keep_rows = [i for i in range(15) if i not in free_rows]
    for ipo in [{}, {"linsol_perturbed_mode": True}]:
      for mi in [1, 3, 6, 500]:
        extra = dict(ipo, max_iter=mi, linsol_iterative_refinement=mi == 500)
        s0 = ca.nlpsol("s0", "ipmc", nlp0, {"structure_detection": "auto", "equality": eq0,
                                            "ipmc": dict(self.IPMC_OPTS, **extra)})
        s1 = ca.nlpsol("s1", "ipmc", nlp1, {"structure_detection": "auto", "equality": eq1,
                                            "ipmc": dict(self.IPMC_OPTS, **extra)})
        r0, r1 = s0(**b0), s1(**b1)
        tag = str(extra)
        self.assertEqual(s0.stats()["iter_count"], s1.stats()["iter_count"], tag)
        self.checkarray(r0["x"], r1["x"], tag+":x", digits=14)
        self.checkarray(r0["lam_x"], r1["lam_x"], tag+":lam_x", digits=14)
        self.checkarray(r0["lam_g"], r1["lam_g"][keep_rows], tag+":lam_g", digits=14)
        for i in free_rows:
          self.assertEqual(float(r1["lam_g"][i]), 0.0, tag)
        if mi == 500:
          self.assertTrue(s1.stats()["success"], tag)

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  def test_ipmc_slacks_lifted_parametric(self):
    self.message("ipmc: lifted slack columns with a parameter")
    # theta goes in f; a parametric penalty is refused
    L1 = lambda s, w: w*ca.sum1(s)
    q = self.slack_case(L1, 1.0, N=20, group_mode="single")
    self.assertEqual(q["ns"], 2)
    theta = ca.SX.sym("theta")
    q["f"] = q["f"] + theta*ca.sumsqr(q["x"])
    q["par"] = theta
    q["nlp"] = dict(q["nlp"], f=q["f"], p=theta)

    seen = []
    for pval in [0.0, 0.05]:
      ref = self.slack_ocp_reference(q, pval=pval)
      for sd in ["none", "manual", "auto"]:
        tag = "theta=%g:%s" % (pval, sd)
        print("test_ipmc_slacks_lifted_parametric", tag)
        solver = self.ipmc_slack_solver(q, sd)
        r = solver(p=pval, **q["bounds"])
        self.assertTrue(solver.stats()["success"], tag)
        # the lifting really is what ran; at "none" nothing is lifted
        self.assertEqual(solver.stats()["n_lift"], 0 if sd == "none" else 2)
        self.check_slack_structure(solver, q, sd)
        for k, d in [("f", 4), ("x", 5), ("g", 6), ("s", 5)]:
          self.checkarray(ref[k], r[k], tag+":"+k, digits=d)
        if sd == "auto":
          seen.append(ca.DM(r["x"]))
          # the parameter through generated C on the lifted path
          if IPMC_CODEGEN:
            self.check_codegen(solver, dict(q["bounds"], p=pval),
                               **IPMC_CODEGEN)
          self.check_serialize(solver, dict(q["bounds"], p=pval))

    # theta moves the solution
    self.assertTrue(float(ca.norm_inf(seen[0]-seen[1])) > 1e-3,
                    "theta did not move the solution")

  @requires_nlpsol("ipmc")
  def test_ipmc_slacks_lifted_distinct_jacobian(self):
    self.message("ipmc: lifted columns on a Jacobian with distinct entries")
    # distinct coefficients expose a permuted nonzero map
    N, nx, nu = 4, 2, 1
    X = [ca.MX.sym("x%d" % k, nx) for k in range(N+1)]
    U = [ca.MX.sym("u%d" % k, nu) for k in range(N)]
    w, g, eq = [], [], []
    f = 0
    for k in range(N):
      w += [X[k], U[k]]
      g.append(X[k+1] - (ca.vertcat(1.1*X[k][0] + 0.3*X[k][1] + 0.7*U[k][0],
                                    0.2*X[k][0] + 0.9*X[k][1])
                         + 0.05*ca.sin(X[k])))
      eq += [True]*nx
      g.append(2*X[k][0] + 3*X[k][1] + 0.1*ca.cos(U[k][0]))
      eq += [False]
      f = f + ca.sumsqr(X[k]) + ca.sumsqr(U[k])
    w.append(X[N])
    f = f + ca.sumsqr(X[N])
    w = ca.vcat(w)
    g = ca.vcat(g)
    ng, nw = g.numel(), w.numel()
    # a lower and an upper budget softening the path row of every stage, so cross-stage
    S = ca.DM(2*(ng+nw), 2)
    for k in range(N):
      S[k*(nx+1)+nx, 0] = 1
      S[ng+nw+k*(nx+1)+nx, 1] = 1
    s = ca.MX.sym("s", 2)
    nlp = {"x": w, "g": g, "s": s, "f_s": 50*ca.sum1(s)+0.5*ca.sumsqr(s), "f": f}
    lbg = ca.DM.zeros(ng)
    ubg = ca.DM.zeros(ng)
    for k in range(N):
      lbg[k*(nx+1)+nx] = -1
      ubg[k*(nx+1)+nx] = 1
    bounds = dict(x0=0.1, lbg=lbg, ubg=ubg, ubs=ca.DM([5, 5]))
    res = {}
    for sd in ["none", "auto"]:
      opts = {"structure_detection": sd, "S": S.sparsity(),
              "print_time": False, "ipmc": dict(self.IPMC_SLACK_OPTS)}
      if sd == "auto":
        opts["equality"] = eq
      solver = ca.nlpsol("solver", "ipmc", nlp, opts)
      r = solver(**bounds)
      self.assertTrue(solver.stats()["success"], sd)
      self.assertEqual(solver.stats()["n_lift"], 0 if sd == "none" else 2)
      res[sd] = r
      if sd == "auto" and IPMC_CODEGEN:
        self.check_codegen(solver, bounds, **IPMC_CODEGEN)
    # at none the columns are stage-local slacks, the reference
    for k, d in [("f", 5), ("x", 5), ("g", 5), ("s", 5)]:
      self.checkarray(res["none"][k], res["auto"][k], "distinct:"+k, digits=d)

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  def test_ipmc_slacks_mixed_columns(self):
    self.message("ipmc: lifted columns and stage-local slacks in one problem")
    L1 = lambda s, w: w*ca.sum1(s)
    p = self.slack_case(L1, 5.0, N=12, group_mode="mixed",
                        soften="upper", v_hi=4.0, p_goal=5.0)
    # one shared pair (stages 1,3,..,11) + one pair per even stage
    self.assertEqual(p["ns"], 14)
    half = p["ns"]//2
    ref = self.slack_ocp_reference(p)
    shared_u = float(ref["s"][half])          # column half is the shared upper one
    local_u = [float(ref["s"][half+j]) for j in range(1, half)]
    n_local = sum(1 for v in local_u if v > 1e-3)
    print("test_ipmc_slacks_mixed_columns shared", shared_u,
          "active per-stage", n_local, "of", half-1)
    self.assertTrue(shared_u > 0.1)           # the lifted column is used
    self.assertTrue(n_local >= 3)             # and so are the stage-local ones

    for sd in ["none", "manual", "auto"]:
      solver = self.ipmc_slack_solver(p, sd)
      r = solver(**p["bounds"])
      self.assertTrue(solver.stats()["success"], sd)
      self.check_slack_structure(solver, p, sd)
      # at none the shared columns are stage-local, so nothing is lifted
      self.assertEqual(solver.stats()["n_lift"], 0 if sd == "none" else 2)
      # ipmc relaxes bounds by 1e-8, ipopt does not: digits 4 to 6
      for k, d in [("f", 4), ("x", 5), ("g", 6), ("s", 4)]:
        self.checkarray(ref[k], r[k], "mixed:"+sd+":"+k, digits=d)
      # helper states and stage-local slacks through serialization together
      self.check_serialize(solver, p["bounds"])

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  def test_ipmc_slacks_lifted_with_nxc(self):
    self.message("ipmc: lifted columns on top of a user-declared nxc")
    q = self.nxc_ocp(N=8, nc=3, p_goal=2.0)
    ns, w = 2, 0.5
    # a lower and an upper budget on the cap of stages 1..N, so lifted
    rows = q["cap_rows"][1:]
    nxu, ngu = q["nlp"]["x"].numel(), q["nlp"]["g"].numel()
    n = ngu+nxu
    S = ca.Sparsity.triplet(2*n, ns, rows+[n+r for r in rows],
                            [0]*len(rows)+[1]*len(rows))
    sv = ca.SX.sym("s", ns)
    p = dict(q, x=q["nlp"]["x"], f=q["nlp"]["f"], g=q["nlp"]["g"],
             nxu=nxu, ngu=ngu, ns=ns, S=S, s=sv, ubs=None,
             nlp={"x": q["nlp"]["x"], "f": q["nlp"]["f"], "g": q["nlp"]["g"],
                  "s": sv, "f_s": w*ca.sum1(sv)})
    ref = self.slack_ocp_reference(p)
    print("test_ipmc_slacks_lifted_with_nxc s", list(numpy.array(ref["s"]).ravel()))
    # non-vacuous: the softened cap really is violated at the optimum
    self.assertTrue(float(ref["s"][1]) > 0.05)

    def solve(sd, nxc):
      opts = {"structure_detection": sd, "print_time": False,
              "ipmc": dict(self.IPMC_SLACK_OPTS), "S": S}
      if sd == "manual":
        opts.update(q["structure"])
      elif sd == "auto":
        opts["equality"] = q["equality"]
      if nxc is not None:
        opts["nxc"] = nxc
      solver = ca.nlpsol("solver", "ipmc", p["nlp"], opts)
      return solver, solver(**q["bounds"])

    # at none the columns are stage-local and nothing is lifted
    sn, rn = solve("none", None)
    self.assertTrue(sn.stats()["success"])
    self.assertEqual(sn.stats()["n_lift"], 0)

    for sd in ["manual", "auto"]:
      s0, r0 = solve(sd, None)          # helpers constant, user's nxc not
      s1, r1 = solve(sd, q["nxc"])      # both trailing blocks constant
      for s_, tag in [(s0, "nxc-off"), (s1, "nxc-on")]:
        self.assertTrue(s_.stats()["success"], sd+":"+tag)
        self.assertEqual(s_.stats()["n_lift"], 2, sd+":"+tag)
      self.assertEqual(s1.stats()["nxc"], q["nxc"])
      self.assertEqual(s0.stats()["nxc"], 0)
      # the same problem, so the same iterate path
      self.assertEqual(s0.stats()["iter_count"], s1.stats()["iter_count"], sd)
      for k, d in self.NXC_DIGITS:
        self.checkarray(r0[k], r1[k], sd+":lift-nxc-vs-lift:"+k, digits=d)
      for k, d in [("f", 4), ("x", 5), ("g", 6), ("s", 4)]:
        self.checkarray(ref[k], r1[k], sd+":lift-nxc-vs-ipopt:"+k, digits=d)
        self.checkarray(rn[k], r1[k], sd+":lift-nxc-vs-none:"+k, digits=d)

    # declaring v constant is still caught with the lifting on
    with self.assertInException("declared constant"):
      solve("manual", q["nc"]+1)

    # the lifted description through codegen and serialization
    sc, _ = solve("auto", q["nxc"])
    if IPMC_CODEGEN:
      self.check_codegen(sc, q["bounds"], **IPMC_CODEGEN)
    self.check_serialize(sc, q["bounds"])

  def linf_chain(self, lifted_by_hand, mode, N=6, n_m=3, nu=2, VB=0.45,
                 w=2.0) -> Dict[str, Any]:
    """Mass chain with an Linf-priced soft speed band, by slack syntax or lifted by hand"""
    nx = 2*n_m + 1
    T = 3.0
    dt = T/N
    xs = ca.SX.sym("xs", nx)
    u = ca.SX.sym("u", nu)
    q, v, p = xs[:n_m], xs[n_m:2*n_m], xs[2*n_m]
    acc = []
    for i in range(n_m):
      qm = q[i-1] if i > 0 else 0
      qp = q[i+1] if i < n_m-1 else 0
      a = (10.0+0.7*i)*(qm-2*q[i]+qp) - (0.2+0.03*i)*v[i] - (0.4+0.05*i)*q[i]**3
      if i < nu: a = a + u[i]
      if i == 0: a = a + 2.0*(p-1.0)
      acc.append(a)
    ode = ca.Function("ode", [xs, u], [ca.vertcat(v, ca.vcat(acc), 0)])
    k1 = ode(xs, u); k2 = ode(xs+dt/2*k1, u); k3 = ode(xs+dt/2*k2, u); k4 = ode(xs+dt*k3, u)
    F = ca.Function("F", [xs, u], [xs + dt/6*(k1+2*k2+2*k3+k4)])
    ncol = n_m if mode == "row" else 1
    col = (lambda i: i) if mode == "row" else (lambda i: 0)
    ne = 2*ncol if lifted_by_hand else 0
    X = [ca.SX.sym("X%d" % k, nx) for k in range(N+1)]
    M = [ca.SX.sym("M%d" % k, ne) for k in range(N+1)]
    U = [ca.SX.sym("U%d" % k, nu) for k in range(N)]
    wv, lbw, ubw, w0, g, lbg, ubg, srow, scol, keep = [], [], [], [], [], [], [], [], [], []
    J = 0.3*(X[0][nx-1]-1.0)**2
    used = [False]*ne
    def var(e, lb, ub, x0, user):
      keep.extend([user]*e.numel())
      wv.append(e); lbw.extend(lb); ubw.extend(ub); w0.extend(x0)
    for k in range(N+1):
      var(X[k], [0]*(nx-1)+[0.8] if k == 0 else [-inf]*nx,
          [0]*(nx-1)+[1.2] if k == 0 else [inf]*nx, [0]*(nx-1)+[1.0], True)
      if ne:
        var(M[k], [0]*ne if k == 0 else [-inf]*ne, [inf]*ne, [0]*ne, False)
      if k < N:
        var(U[k], [-25]*nu, [25]*nu, [0]*nu, True)
        g.append(X[k+1]-F(X[k], U[k])); lbg += [0]*nx; ubg += [0]*nx
        if ne:
          g.append(M[k+1]-M[k]); lbg += [0]*ne; ubg += [0]*ne
        J = J + dt*ca.sumsqr(U[k])
      if k >= 1:
        for i in range(n_m):
          c = 1.3*v[i] + 0.4*q[i] + 0.2*q[i]**2
          c = ca.substitute(c, xs, X[k])
          lo = -inf if i == 0 else -VB
          if not ne:
            srow.append(len(lbg)); scol.append(col(i))
            g.append(c); lbg.append(lo); ubg.append(VB)
            continue
          if lo > -inf:
            e = 2*col(i); used[e] = True
            g.append(c+M[k][e]); lbg.append(lo); ubg.append(inf)
          e = 2*col(i)+1; used[e] = True
          g.append(c-M[k][e]); lbg.append(-inf); ubg.append(VB)
      if k == N:
        g.append(X[N][0]); lbg.append(0.6); ubg.append(inf)
    x = ca.vcat(wv)
    gg = ca.vcat(g)
    eq = [a == b for a, b in zip(lbg, ubg)]
    bounds = dict(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
    if ne:
      nlp = {"x": x, "g": gg, "f": J + w*ca.sum1(M[0])}
      return dict(nlp=nlp, bounds=bounds, eq=eq, keep=numpy.array(keep), ne=ne, nx=nx)
    # a lower budget (column 2c) and an upper budget (2c+1) per column of the twin
    s = ca.SX.sym("s", 2*ncol)
    n = gg.numel()+x.numel()
    S = ca.Sparsity.triplet(2*n, 2*ncol, srow+[n+r for r in srow],
                            [2*c for c in scol]+[2*c+1 for c in scol])
    nlp = {"x": x, "g": gg, "f": J, "s": s, "f_s": w*ca.sum1(s)}
    return dict(nlp=nlp, bounds=bounds, eq=eq, S=S, ne=2*ncol, nx=nx)

  @requires_nlpsol("ipmc")
  def test_ipmc_slacks_lifted_linf_twin(self):
    self.message("ipmc: lifted Linf slacks against the hand-lifted twin, step for step")
    for mode in ["row", "side"]:
      hand = self.linf_chain(True, mode)
      nat = self.linf_chain(False, mode)
      for ipo in [{}, {"linsol_perturbed_mode": True}]:
        def solve(prob, extra, native):
          opts = {"structure_detection": "auto", "equality": prob["eq"], "print_time": False,
                  "ipmc": dict(self.IPMC_SLACK_OPTS, **dict(ipo, **extra))}
          if native:
            opts.update(S=prob["S"], nxc=1)
          solver = ca.nlpsol("solver", "ipmc", prob["nlp"], opts)
          return solver, solver(**prob["bounds"])
        tag = mode + str(ipo)
        sh, rh = solve(hand, {}, False)
        sn, rn = solve(nat, {}, True)
        print("test_ipmc_slacks_lifted_linf_twin", tag, sh.stats()["iter_count"], float(rn["f"]))
        self.assertTrue(sh.stats()["success"] and sn.stats()["success"], tag)
        # one helper per column of S, the user's nxc on top
        self.assertEqual(sn.stats()["n_lift"], nat["ne"], tag)
        self.assertEqual(sn.stats()["nxc"], 1, tag)
        self.assertEqual(sh.stats()["iter_count"], sn.stats()["iter_count"], tag)
        self.checkarray(rh["f"], rn["f"], tag+":f", digits=12)
        xh = numpy.array(rh["x"]).ravel()[hand["keep"]]
        self.checkarray(ca.DM(xh), rn["x"], tag+":x", digits=12)
        # the budgets are the stage-0 helper states of the twin
        mh = numpy.array(rh["x"]).ravel()[~hand["keep"]][:hand["ne"]]
        self.checkarray(ca.DM(mh), rn["s"], tag+":s", digits=10)
        # non-vacuous: the band is violated, so the budgets are live
        self.assertTrue(float(ca.norm_inf(rn["s"])) > 1e-2, tag)
        # truncated solves without refinement expose a wrong factorization
        for mi in [1, 3, 6]:
          extra = {"max_iter": mi, "linsol_iterative_refinement": False}
          _, th = solve(hand, extra, False)
          _, tn = solve(nat, extra, True)
          dx = float(numpy.max(numpy.abs(numpy.array(th["x"]).ravel()[hand["keep"]]
                                         - numpy.array(tn["x"]).ravel())))
          self.assertTrue(dx < 1e-9, "%s max_iter %d: dx %g" % (tag, mi, dx))
      if IPMC_CODEGEN:
        ipo = {}
        sc, _ = solve(nat, {}, True)
        self.check_codegen(sc, nat["bounds"], **IPMC_CODEGEN)
      self.check_serialize(sn, nat["bounds"])

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  @memory_heavy()
  def test_ipmc_slacks_ubs_finite(self):
    self.message("ipmc: slacks saturating at a finite ubs")
    L1 = lambda s, w: w*ca.sum1(s)
    cases = [
      # name, w, kwargs, (ubsl, ubsu), min saturated (s_l, s_u), hard twin feasible
      # a two-sided simple bound (S_x): both halves saturate
      ("ubs_bound_band_x_N10", 0.5,
       dict(N=10, soften="bound_band_x", v_lo=3.0, v_hi=6.0),
       (1.0, 0.6), (1, 2), True),
      # softened equality rows: eleven of twenty upper slacks saturate
      ("ubs_equality_N20", 0.5, dict(N=20, v_lo=3.0, v_hi=3.0),
       (1.8, 3.0), (1, 11), False),
      # lifted columns: ubs bounds the helper states
      ("ubs_linf_N20", 1.0, dict(N=20, group_mode="single"),
       (0.1, 0.1), (1, 1), True),
    ]
    for name, w, kwargs, cap, sat_min, hard_feasible in cases:
      ns = 1 if kwargs.get("group_mode") == "single" else kwargs["N"]
      ubs = ca.vertcat(cap[0]*ca.DM.ones(ns), cap[1]*ca.DM.ones(ns))
      p = self.slack_case(L1, w, ubs=ubs, **kwargs)
      self.assertEqual(p["ns"], 2*ns)
      ref = self.slack_ocp_reference(p)
      rh, hard_ok = self.ipmc_hard_solve(p)
      self.assertEqual(hard_ok, hard_feasible)

      def saturated(s):
        return (sum(1 for i in range(ns) if float(s[i]) > cap[0]-1e-6),
                sum(1 for i in range(ns) if float(s[ns+i]) > cap[1]-1e-6))
      sat = saturated(ref["s"])
      print("test_ipmc_slacks_ubs_finite", name, "saturated", sat,
            "of", ns, "f", float(ref["f"]))
      self.assertTrue(sat[0] >= sat_min[0] and sat[1] >= sat_min[1])
      # ipmc relaxes bounds by 1e-8: a saturated slack sits ~3e-8 above ubs
      self.assertTrue(float(ca.mmax(ref["s"]-ubs)) < 1e-6)

      for sd in ["none", "manual", "auto"]:
        print("test_ipmc_slacks_ubs_finite", name, sd)
        solver = self.ipmc_slack_solver(p, sd)
        r = solver(ubs=ubs, **p["bounds"])
        self.assertTrue(solver.stats()["success"])
        self.check_slack_structure(solver, p, sd)
        # the native path saturates the same slacks, and never overshoots
        self.assertEqual(saturated(r["s"]), sat)
        self.assertTrue(float(ca.mmax(r["s"]-ubs)) < 1e-6)
        # most active bounds in the file: g digits 5
        for k, d in [("f", 4), ("x", 5), ("g", 5), ("s", 5)]:
          self.checkarray(ref[k], r[k], name+":"+sd+":"+k, digits=d)
        # ubs is a runtime input, fed to the round trip explicitly
        self.check_serialize(solver, dict(p["bounds"], ubs=ubs))

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  def test_ipmc_slacks_lifted_ubs(self):
    self.message("ipmc: ubs on lifted columns, changing between solves; 0 hardens the rows")
    L1 = lambda s, w: w*ca.sum1(s)
    p = self.slack_case(L1, 1.0, N=20, group_mode="single")
    self.assertEqual(p["ns"], 2)
    rh, hard_ok = self.ipmc_hard_solve(p)
    self.assertTrue(hard_ok)
    calls = [("free", [inf, inf]), ("lower capped", [0.05, inf]),
             ("upper capped", [inf, 0.05]), ("both zero", [0, 0]),
             ("lower zero", [0, inf]), ("upper zero", [inf, 0]),
             ("both capped", [0.05, 0.05]), ("both zero again", [0, 0])]
    # at none the columns are stage-local, where ipmc needs a positive bound
    with self.assertInException("has the upper bound ubs = 0, but relaxes a finite bound"):
      self.ipmc_slack_solver(p, "none")(ubs=[0, inf], **p["bounds"])
    for sd in ["manual", "auto"]:
      reused = self.ipmc_slack_solver(p, sd)
      for name, ubs in calls:
        print("test_ipmc_slacks_lifted_ubs", name, sd)
        r = reused(ubs=ubs, **p["bounds"])
        self.assertTrue(reused.stats()["success"])
        self.assertEqual(reused.stats()["n_lift"], 2)
        self.check_slack_structure(reused, p, sd)
        # a reused solver gives bit for bit what a fresh one gives
        fresh = self.ipmc_slack_solver(p, sd)
        rf = fresh(ubs=ubs, **p["bounds"])
        for k in ["f", "x", "s", "g", "lam_x", "lam_g", "lam_s"]:
          self.assertTrue(numpy.array_equal(numpy.array(r[k]), numpy.array(rf[k])),
                          name+":"+sd+":"+k+": reused != fresh")
        if ubs == [0, 0]:
          # nothing can be relaxed: the hard problem
          self.checkarray(ca.DM.zeros(2), r["s"], name+":"+sd+":s", digits=15)
          for k in ["f", "x", "g"]:
            self.checkarray(rh[k], r[k], name+":"+sd+":"+k, digits=7)
          continue
        ref = self.slack_ocp_reference(dict(p, ubs=ubs))
        if name == "free":
          self.assertTrue(float(ref["s"][0]) > 0.1 and float(ref["s"][1]) > 0.1)
        for i in range(2):
          if ubs[i] == 0:
            self.assertEqual(float(r["s"][i]), 0.0)
          elif ubs[i] < inf:
            # the cap binds, and its multiplier is reported
            self.assertTrue(abs(float(r["s"][i]) - ubs[i]) < 1e-6)
            self.assertTrue(float(r["lam_s"][i]) > 1e-3)
        for k, d in [("f", 4), ("x", 5), ("g", 5), ("s", 5), ("lam_s", 4)]:
          self.checkarray(ref[k], r[k], name+":"+sd+":"+k, digits=d)
      self.check_serialize(reused, dict(p["bounds"], ubs=[0, 0]))
      if IPMC_CODEGEN:
        self.check_codegen(reused, dict(p["bounds"], ubs=[0, inf]), **IPMC_CODEGEN)

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  def test_ipmc_slacks_equality_row(self):
    self.message("ipmc: a softened equality row")
    # v is pinned to 3 but the goal needs a mean speed of 5, so the hard twin fails
    p = self.slack_case(lambda s, w: w*ca.sum1(s), 0.5, N=20,
                        v_lo=3.0, v_hi=3.0)
    # every softened row is a declared equality row
    for r in set(i % (p["ngu"]+p["nxu"]) for i in p["S"].row()):
      self.assertTrue(r < p["ngu"])
      self.assertEqual(p["lbg"][r], p["ubg"][r])
      self.assertTrue(p["equality"][r])

    ref = self.slack_ocp_reference(p)
    rh, hard_ok = self.ipmc_hard_solve(p)
    self.assertFalse(hard_ok)
    nsl, nsu = self.slack_activity(ref["s"], p)
    print("test_ipmc_slacks_equality_row active", (nsl, nsu),
          "f", float(ref["f"]))
    self.assertTrue(nsl >= 4 and nsu >= 16)

    for sd in ["none", "manual", "auto"]:
      print("test_ipmc_slacks_equality_row", sd)
      solver = self.ipmc_slack_solver(p, sd)
      r = solver(**p["bounds"])
      self.assertTrue(solver.stats()["success"])
      self.check_slack_structure(solver, p, sd)
      for k, d in [("f", 4), ("x", 5), ("g", 6), ("s", 5)]:
        self.checkarray(ref[k], r[k], "equality_row_N20:"+sd+":"+k, digits=d)

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  def test_ipmc_slacks_asymmetric_penalty(self):
    self.message("ipmc: different weights on s_l and s_u")
    def asym(wl, wu):
      return lambda s, w: wl*ca.sum1(s[:s.numel()//2]) \
                          + wu*ca.sum1(s[s.numel()//2:])
    kwargs = dict(N=20, soften="bound_band_x", v_lo=3.0, v_hi=6.0)
    sums = {}
    for wl, wu in [(0.3, 3.0), (3.0, 0.3)]:
      name = "asym_%g_%g" % (wl, wu)
      p = self.slack_case(asym(wl, wu), None, **kwargs)
      ref = self.slack_ocp_reference(p)
      rh, hard_ok = self.ipmc_hard_solve(p)
      self.assertFalse(hard_ok)      # non-vacuity: the hard twin fails
      ns = p["ns"]//2
      nsl, nsu = self.slack_activity(ref["s"], p)
      print("test_ipmc_slacks_asymmetric_penalty", name,
            "active", (nsl, nsu), "f", float(ref["f"]))
      # both halves have to be active or the asymmetry is untested
      self.assertTrue(nsl >= 3 and nsu >= 9)
      sums[(wl, wu)] = (float(ca.sum1(ref["s"][:ns])),
                        float(ca.sum1(ref["s"][ns:])), float(ref["f"]))

      for sd in ["none", "manual", "auto"]:
        print("test_ipmc_slacks_asymmetric_penalty", name, sd)
        solver = self.ipmc_slack_solver(p, sd)
        r = solver(**p["bounds"])
        self.assertTrue(solver.stats()["success"])
        self.check_slack_structure(solver, p, sd)
        for k, d in [("f", 4), ("x", 5), ("g", 6), ("s", 5)]:
          self.checkarray(ref[k], r[k], name+":"+sd+":"+k, digits=d)

    # the run with the cheap upper slack uses more of it
    self.assertTrue(sums[(3.0, 0.3)][1] > sums[(0.3, 3.0)][1] + 1.0)
    self.assertTrue(abs(sums[(3.0, 0.3)][2] - sums[(0.3, 3.0)][2]) > 1.0)

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  def test_ipmc_slacks_parametric_penalty(self):
    self.message("ipmc: a slack penalty weighted by p, retuned between solves")
    par = ca.SX.sym("w_slack")
    p = self.slack_case(lambda s, w: par*ca.sum1(s), None, par=par,
                        N=20, soften="bound_band_x", v_lo=3.0, v_hi=6.0)
    weights = [0.3, 3.0]
    ref = {pv: self.slack_ocp_reference(p, pv) for pv in weights}
    tot = {pv: float(ca.sum1(ref[pv]["s"])) for pv in weights}
    print("test_ipmc_slacks_parametric_penalty sum|s|", tot)
    # the two weights give different answers, the cheap one relaxes more
    self.assertTrue(tot[0.3] > tot[3.0] + 1.0)
    self.assertTrue(float(ca.norm_inf(ref[0.3]["x"]-ref[3.0]["x"])) > 0.1)
    for pv in weights:
      nsl, nsu = self.slack_activity(ref[pv]["s"], p)
      self.assertTrue(nsl >= 1 and nsu >= 9)

    for sd in ["none", "manual", "auto"]:
      print("test_ipmc_slacks_parametric_penalty", sd)
      solver = self.ipmc_slack_solver(p, sd)
      # one solver object; the third solve repeats the first, which a stale weight cannot
      out = []
      for pv in [0.3, 3.0, 0.3]:
        r = solver(p=pv, **p["bounds"])
        self.assertTrue(solver.stats()["success"])
        self.check_slack_structure(solver, p, sd)
        for k, d in [("f", 4), ("x", 5), ("g", 6), ("s", 5)]:
          self.checkarray(ref[pv][k], r[k], "%g:%s:%s" % (pv, sd, k), digits=d)
        out.append(r)
      for k in ["f", "x", "g", "s"]:
        self.checkarray(out[0][k], out[2][k], "repeat:"+sd+":"+k, digits=12)
      if sd == "auto":
        if IPMC_CODEGEN:
          self.check_codegen(solver, dict(p["bounds"], p=3.0), **IPMC_CODEGEN)
        self.check_serialize(solver, dict(p["bounds"], p=3.0))

    # expand_slacks tracks p
    opts = {"structure_detection": "none", "ipmc": dict(self.IPMC_SLACK_OPTS), "S": p["S"],
            "expand_slacks": True}
    solver = ca.nlpsol("solver", "ipmc", p["nlp"], opts)
    out = []
    for pv in [0.3, 3.0, 0.3]:
      r = solver(p=pv, **p["bounds"])
      self.assertTrue(solver.stats()["success"])
      for k, d in [("f", 4), ("x", 5), ("g", 6), ("s", 5)]:
        self.checkarray(ref[pv][k], r[k],
                        "parametric_%g:expand:%s" % (pv, k), digits=d)
      out.append(r)
    # the third solve reproduces the first
    for k in ["f", "x", "g", "s"]:
      self.checkarray(out[0][k], out[2][k],
                      "parametric_repeat:expand:"+k, digits=12)

  @requires_nlpsol("ipmc")
  def test_ipmc_no_convexify(self):
    self.message("ipmc refuses convexify_strategy, which it never applied")
    nlp, b, eq = self.ipmc_small_ocp()
    for o in ["convexify_strategy", "convexify_margin"]:
      with self.assertInException(o):
        ca.nlpsol("s", "ipmc", nlp, {o: {"convexify_strategy": "regularize",
                                         "convexify_margin": 1e-7}[o]})

  @requires_nlpsol("ipmc")
  def test_ipmc_return_status(self):
    self.message("ipmc reports its status by name")
    nlp, b, eq = self.ipmc_small_ocp()
    for mi, status in [(500, "IPMC_SOLVED"), (1, "IPMC_FAILED")]:
      solver = ca.nlpsol("s", "ipmc", nlp, {"structure_detection": "auto", "equality": eq,
                                           "ipmc": dict(self.IPMC_OPTS, max_iter=mi)})
      solver(**b)
      self.assertEqual(solver.stats()["return_status"], status)

  @requires_nlpsol("ipmc")
  def test_ipmc_codegen_two_solvers(self):
    self.message("two ipmc solvers in one generated file")
    if not IPMC_CODEGEN: return
    nlp, b, eq = self.ipmc_small_ocp()
    sols = [ca.nlpsol(n, "ipmc", nlp, dict({"structure_detection": sd,
                                            "ipmc": dict(self.IPMC_OPTS)}, **extra))
            for n, sd, extra in [("ipmc_two_a", "auto", {"equality": eq}),
                                 ("ipmc_two_b", "none", {})]]
    cg = ca.CodeGenerator("ipmc_two.c")
    for s in sols: cg.add(s)
    cg.generate()
    r = self.compile_external("ipmc_two_a", "ipmc_two.c", **IPMC_CODEGEN)
    if r is not None:
      for s in sols:
        F = ca.external(s.name(), r[1])
        ref = s(**b)
        out = F(**b)
        for k in ["x", "f", "lam_x", "lam_g"]:
          self.checkarray(ref[k], out[k], s.name()+":"+k, digits=12)
      os.remove(r[1])
    os.remove("ipmc_two.c")

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  def test_ipmc_slacks_one_sided(self):
    self.message("ipmc: ubs = 0 on a slack that relaxes no finite bound is inert")
    nlp, b, eq = self.ipmc_small_ocp()
    S = ca.Sparsity.triplet(2*(6+7), 6, [1, 3, 5, 13+1, 13+3, 13+5], range(6))
    s = ca.SX.sym("s", 6)
    solver = ca.nlpsol("s", "ipmc", dict(nlp, s=s, f_s=0.1*ca.sum1(s)),
                       {"structure_detection": "auto", "equality": eq, "S": S,
                        "print_time": False, "ipmc": dict(self.IPMC_OPTS)})
    # the softened rows have no lower bound; the upper one binds
    args = dict(b, ubg=[0, -0.5, 0, -0.5, 0, -0.5])
    ref = solver(**args)
    r = solver(**dict(args, ubs=[0, 0, 0, inf, inf, inf]))
    self.assertTrue(solver.stats()["success"])
    self.checkarray(ca.DM.zeros(3), r["s"][:3], "s_l", digits=15)
    self.checkarray(ca.DM.zeros(3), r["lam_s"][:3], "lam_s_l", digits=15)
    self.assertTrue(float(ca.norm_inf(r["s"][3:])) > 0.1)
    for k in ["x", "f", "lam_g", "lam_x"]:
      self.checkarray(ref[k], r[k], k, digits=6)
    # a side with a finite bound is still refused
    with self.assertInException("slack column 4 has the upper bound ubs = 0, but relaxes a finite bound"):
      solver(**dict(args, ubs=[inf, inf, inf, inf, 0, inf]))
    if IPMC_CODEGEN:
      self.check_codegen(solver, dict(args, ubs=[0, 0, 0, inf, inf, inf]), **IPMC_CODEGEN)
    # Opti makes one column per slack element, on the side the user wrote
    opti = ca.Opti()
    X = opti.variable(4); U = opti.variable(3); v = opti.slack(3)
    opti.subject_to(X[0] == 0)
    for i in range(3):
      opti.subject_to(X[i+1] == X[i] + U[i])
      opti.subject_to(opti.bounded(-1, U[i], 1))
      opti.subject_to(X[i+1] <= 0.5 + v[i])
    opti.minimize(ca.sumsqr(X-3) + ca.sumsqr(U) + 2*ca.sum1(v))
    f = {}
    for name, o in [("ipopt", {"ipopt.print_level": 0, "ipopt.tol": 1e-10}),
                    ("ipmc", {"ipmc": dict(self.IPMC_OPTS)})]:
      opti.solver(name, dict(o, print_time=False))
      f[name] = opti.solve().value(opti.f)
    self.checkarray(f["ipopt"], f["ipmc"], "opti f", digits=6)

if __name__ == '__main__':
    unittest.main()
