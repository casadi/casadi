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
from helpers import *

import os

# Codegen settings for the ipmc entry, mirroring nlp.py's: ipmc is not
# installed into casadi's tree, so the generated code has to be pointed at
# its build directory, and the check disables itself -- with a printed
# reason -- when either that or casadi's own include tree is missing.  The
# latter is absent unless casadi was configured WITH_SELFCONTAINED, which is
# the same thing that silently disables the fatrop codegen check below.
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
        
  
  # =====================================================================
  # ipmc on a real staircase OCP.
  #
  # nlp.py runs ipmc over the whole solver battery, but always at the
  # default structure_detection="none", where the entire problem is one
  # dense stage and the Riccati recursion has nothing to recurse over.
  # Everything this plugin exists for -- the block detection, the
  # A B I / C D projections, the stagewise equality/inequality split and
  # the readback of a stagewise multiplier vector into casadi's flat one --
  # is only reached once a horizon is actually detected.
  #
  # Ground truth is ipopt on the very same NLP, with bound_relax_factor=0
  # so that the reference is the optimum of the problem as written: ipmc,
  # like ipopt at its defaults, widens every bound by
  # bound_relax_factor*max(1,|bound|) = 1e-8 before starting the barrier
  # and never undoes it, so the point it returns sits ~7e-8 outside the
  # feasible set, and its objective is correspondingly a little BELOW the
  # true optimum.  That is what sets the tolerances below: measured worst
  # case over the three settings is f 7.6e-6, x 1.9e-6, g 2.4e-7 and
  # 1.5e-5 on both multiplier vectors, and the digits are the tightest
  # value leaving >= ~2x headroom on that.
  # =====================================================================
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
    """A minimum-effort sprint OCP, laid out exactly as ipmc wants.

        p' = v,   v' = u - 0.02 v^2        RK4, N intervals, horizon T

    Decision vector   x = [x_0, u_0, x_1, u_1, ..., x_{N-1}, u_{N-1}, x_N]
    Constraint vector g = [dyn_0, path_0, ..., dyn_{N-1}, path_{N-1}, path_N]

    so nx = [2]*(N+1), nu = [1]*N + [0], ng = [1]*N + [2] and the Jacobian is
    the staircase.  The drag term keeps the Hessian non-constant, so the
    recursion really is re-run at every iteration; the terminal row
    p_N >= p_goal is what forces the speed limit to bind.
    """
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
    """The three structure_detection settings on a staircase OCP.

    "none" is the dense fallback and is what nlp.py already covers; it is
    kept here as the control.  "manual" takes the partition from the user
    and "auto" has to find it in the constraint-Jacobian sparsity, and the
    point of the test is that all three land on the same answer as a dense
    ipopt solve of the very same NLP -- primal AND dual, since the
    multipliers are the part that has to be scattered back out of ipmc's
    stagewise layout.

    stats() is asserted as well: a mis-detected horizon that still happened
    to converge would otherwise pass silently.
    """
    p = self.sprint_ocp()
    ref = ca.nlpsol("reference", "ipopt", p["nlp"],
                    dict(self.IPOPT_HARD_REF_OPTS,
                         ipopt=dict(self.IPOPT_HARD_REF_OPTS["ipopt"])))
    r_ref = ref(**p["bounds"])
    self.assertTrue(ref.stats()["success"])
    # Non-vacuity: the speed cap really binds, so there is an active
    # inequality on almost every stage and the dual comparison below has
    # something to compare.
    v = [float(r_ref["x"][3*k+1]) for k in range(1, p["structure"]["N"]+1)]
    n_active = sum(1 for vk in v if vk > 6.0-1e-6)
    print("test_ipmc_structured active speed rows", n_active, "f", float(r_ref["f"]))
    self.assertTrue(n_active >= 5)

    for sd in ["none", "manual", "auto"]:
      print("test_ipmc_structured", sd)
      opts = {"structure_detection": sd, "equality": p["equality"],
              "ipmc": dict(self.IPMC_OPTS)}
      if sd == "manual":
        opts.update(p["structure"])
      solver = ca.nlpsol("solver", "ipmc", p["nlp"], opts)
      r = solver(**p["bounds"])
      self.assertTrue(solver.stats()["success"])

      if sd != "none":
        # the horizon that was used is the one the user's problem has
        st = solver.stats()
        self.assertEqual(st["N"], p["structure"]["N"])
        for k in ["nx", "nu", "ng"]:
          self.checkarray(ca.DM(p["structure"][k]), ca.DM(st[k]),
                          "structure:"+k, digits=12)

      for k, d in [("f", 4), ("x", 5), ("g", 6), ("lam_g", 4), ("lam_x", 4)]:
        self.checkarray(r_ref[k], r[k], sd+":"+k, digits=d)

    # One structured solver through codegen and serialization: the block
    # partition and the offsets have to survive both.
    solver = ca.nlpsol("solver", "ipmc", p["nlp"],
                       dict({"structure_detection": "auto",
                             "equality": p["equality"],
                             "ipmc": dict(self.IPMC_OPTS)}))
    if IPMC_CODEGEN:
      self.check_codegen(solver, p["bounds"], **IPMC_CODEGEN)
    self.check_serialize(solver, p["bounds"])

  # =====================================================================
  # nxc -- constant states
  #
  # A constant state is one whose dynamics are the identity, x_{k+1}=x_k:
  # a parameter that the solver is allowed to choose, carried along the
  # horizon.  Written out in full it costs a state dimension everywhere,
  # and the Riccati recursion is cubic in that dimension.  The nxc option
  # says "the last nxc states of every stage are of that kind", the
  # interface then hands ipmc a dynamics block with those rows cut off,
  # and ipmc substitutes the identity symbolically instead of
  # multiplying it out.  The problem the user writes does not change at
  # all -- that is what the tests below check.
  # =====================================================================

  # See test_ipmc_nxc's docstring for where these come from.
  NXC_DIGITS = [("f", 12), ("x", 12), ("g", 13), ("lam_g", 10), ("lam_x", 10)]
  # ... and test_ipmc_nxc_restoration's, for the ones that hold on a solve
  # that ends in the restoration phase rather than at a solution.
  NXC_RESTO_DIGITS = [("f", 6), ("x", 8), ("g", 8), ("lam_g", 6), ("lam_x", 6)]

  def nxc_ocp(self, N=10, nc=4, T=1.0, p_goal=2.0, u_max=40.0, W=10.0,
              bound_every=True, terminal_eq=False) -> Dict[str, Any]:
    """A sprint OCP whose state carries nc CONSTANT components.

        x_k = [p_k, v_k, c_0 .. c_{nc-1}]        nx = 2+nc at every stage
        p' = v,   v' = u - c_0 v^2               <- c_0 is IN THE DYNAMICS
        v_k <= c_1 + sum_{j>=2} c_j              <- every stage, so the
                                                    constants sit in a path
                                                    constraint that BINDS
        p_N >= p_goal                            <- what forces the effort
        0.005 <= c_0 <= 0.05, 0.5 <= c_1 <= 20, 0 <= c_j <= 0.4

    objective  dt*sum u^2 + W c_1^2 + sum_j W_j c_j^2, so raising the cap
    buys a cheaper trajectory and the constants are priced against it.

    Everything the constants do here is deliberate.  They enter the
    dynamics (c_0), the objective (all), a path constraint at every stage
    (c_1 and up) and their own bounds, and at the solution c_0 sits on its
    lower bound, at least one c_j on its upper bound and the cap is active
    on most stages -- so a test that compares multipliers has multipliers
    to compare.  A formulation where the constants were inert would pass
    every check below while proving nothing.

    The dynamics rows are written as the extended form the option is
    about:  x_{k+1}[:2] = F(x_k, u_k)  and  x_{k+1}[2:] = x_k[2:], the
    second block being the [0 I] that nxc lets the interface cut away.
    """
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
      # A bound on a constant only has to be imposed once.  Imposing it on
      # every stage is what a user would write and is covered by default;
      # bound_every=False is the same problem with the redundant copies
      # dropped, which is how the timings avoid being dominated by the
      # inequality rows they create.
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
    # row index of the v_k <= cap row of every stage, for the tests that
    # soften it (test_ipmc_slacks_lifted_with_nxc)
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
      # Two EQUALITY rows at the terminal stage, both involving the constants.
      # The terminal stage has no input to pivot on, so neither row can be
      # eliminated there: they are pushed back to stage N-1 as the H block of
      # Eq. (3.9), which is the only way to reach the second structured path --
      # the one that pulls H through the dynamics.  Without them gamma_k is 0
      # at every stage and that code never runs.
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
    opts = {"structure_detection": sd, "equality": p["equality"],
            "print_time": False, "ipmc": dict(self.IPMC_OPTS)}
    if sd == "manual":
      opts.update(p["structure"])
    if nxc is not None:
      opts["nxc"] = nxc
    if opts_extra:
      opts.update(opts_extra)
    solver = ca.nlpsol("solver", "ipmc", p["nlp"], opts)
    return solver, solver(**p["bounds"])

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  def test_ipmc_nxc(self):
    """nxc must not change the answer -- primal, dual or objective.

    Three comparisons, on the same OCP, at both structure_detection
    settings the option applies to:

      1. nxc absent versus nxc all-zero.  BIT identical is asserted, not
         approximate: an all-zero declaration must take exactly the code
         path that was there before, and anything else means the "absent
         is free" claim is not true.
      2. nxc absent versus nxc = the real thing.  This is the claim.  The
         two solves are not merely close at the solution, they take the
         SAME NUMBER OF ITERATIONS, which is a far sharper test: the
         factorization feeds the line search, so a wrong Riccati
         recursion moves the iterate path long before it moves the
         answer.
      3. both against a dense ipopt solve, so that "the same" cannot mean
         "the same and wrong".

    TOLERANCES.  Measured over the 37 converging configurations these
    tests cover -- N in 1..10, nc in 2..4, every nxc in 0..nc, both
    structure_detection settings, both terminal-constraint shapes and
    the perturbed (delta_c) factorization -- the worst nxc-on vs nxc-off
    difference in checkarray's own metric was |df| 3.6e-14, |dx| 4.1e-14,
    |dg| 2.0e-15, |dlam_g| 2.0e-11, |dlam_x| 1.9e-11, and the iteration
    count was equal in every single one.  checkarray compares against
    0.5e-digits (absolute below 1e3, on the leading digit above it), so
    the values below are the tightest leaving >= 2x headroom: f 12 (14x),
    x 12 (12x), g 13 (25x), the multipliers 10 (2.4x).  The multipliers
    are the loose ones because they are what the factorization produces
    most directly.  Against ipopt the tolerances are the looser ones
    test_ipmc_structured already documents, for the same reason
    (bound_relax_factor).

    test_ipmc_nxc_restoration measures its own, much looser, set: see
    NXC_RESTO_DIGITS.
    """
    p = self.nxc_ocp(N=10, nc=4)
    N, nc, nx = p["N"], p["nc"], 2+p["nc"]

    ref = ca.nlpsol("reference", "ipopt", p["nlp"],
                    dict(self.IPOPT_HARD_REF_OPTS,
                         ipopt=dict(self.IPOPT_HARD_REF_OPTS["ipopt"])))
    r_ref = ref(**p["bounds"])
    self.assertTrue(ref.stats()["success"])

    # ---- non-vacuity: the constants have to be doing something ----
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
    # at least one of the others is pinned on its upper bound, and at least
    # one is strictly inside it -- so both kinds of constant are covered
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

      # 1. all-zero nxc is the old code path, exactly
      for k in ["f", "x", "g", "lam_g", "lam_x"]:
        self.assertEqual(float(ca.norm_inf(r0[k]-rz[k])), 0.0,
                         sd+": nxc=0 is not bit-identical on "+k)

      # 2. the same iterate path, and the same answer
      self.assertEqual(s0.stats()["iter_count"], s1.stats()["iter_count"])
      for k, d in self.NXC_DIGITS:
        self.checkarray(r0[k], r1[k], sd+":nxc-vs-plain:"+k, digits=d)

      # 3. and it is the right answer
      for k, d in [("f", 4), ("x", 5), ("g", 6), ("lam_g", 4), ("lam_x", 4)]:
        self.checkarray(r_ref[k], r1[k], sd+":nxc-vs-ipopt:"+k, digits=d)

  @requires_nlpsol("ipmc")
  def test_ipmc_nxc_horizons(self):
    """Short horizons, and declaring only part of the constant block.

    N = 1 is the reason this test exists: it is the one horizon with a
    single dynamics transition and no interior stage, so the backward
    recursion is only its first and last step and every offset in the
    substitution has to be right on a step that is simultaneously both.
    N = 2 adds exactly one interior stage.

    The nxc sweep is the second axis.  nxc need not name every constant
    the problem has -- declaring fewer is legal, exploits less, and must
    still land on the same answer, which is what walks the trailing
    offset crow = nu+nx-nxc over values other than "all of them".  nxc =
    0 is in the sweep because it is the promise that the option is free
    when it is not used, and it is asserted BIT-identical, not close.
    """
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
          for k, d in self.NXC_DIGITS:
            self.checkarray(r0[k], r1[k], tag+":"+k, digits=d)

  @requires_nlpsol("ipmc")
  def test_ipmc_nxc_paths(self):
    """Both structured branches, and the delta_c factorization.

    The dynamics substitution has two structured branches, not one, and
    an ordinary OCP only ever reaches the first.  The second is the one
    that pulls stage k+1's leftover equality rows back through the
    dynamics (the H block of Eq. (3.9)): it runs only when a stage could
    not eliminate all its equality rows on its own inputs, and gamma_k is
    identically zero -- so that code never executes at all -- unless the
    problem has stagewise equalities that outnumber the inputs.
    terminal_eq=True is what creates them: two equality rows at the
    terminal stage, which has no input to pivot on, so both are pushed
    back.  This is not a hypothetical gap.  Writing that branch from a
    comment that had the transposed-add's m and n the wrong way round
    produced a factorization that was wrong from the very first dual
    estimate, converged anyway, and would have gone unnoticed by every
    other test here.

    linsol_perturbed_mode is the second axis: it routes every solve
    through the degenerate (delta_c) factorization instead of the normal
    one.  That is a SECOND copy of the substitution, with its own call to
    the shared helper, and nothing else in this file reaches it.
    """
    for N, nc in [(5, 3), (10, 4)]:
      for te in [False, True]:
        p = self.nxc_ocp(N=N, nc=nc, terminal_eq=te)
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
    """nxc is a declaration, so every way of declaring nonsense must fail.

    The last one is the one that matters.  nxc is not detected, it is
    promised, and a promise that the dynamics do not keep would otherwise
    hand ipmc a Jacobian block that is not the problem's -- a wrong
    Newton step, converging to something, silently.  The interface checks
    the rows it drops against the identity once per solve and refuses.
    """
    p = self.nxc_ocp(N=6, nc=3)
    N, nc = p["N"], p["nc"]

    # structure_detection none has no horizon for nxc to describe -- and
    # even nxc = 0, the value that changes nothing, is refused there,
    # because declaring it at all says the user thinks there is a horizon
    with self.assertInException("structure_detection"):
      ca.nlpsol("solver", "ipmc", p["nlp"],
                {"structure_detection": "none", "equality": p["equality"],
                 "nxc": 0, "ipmc": dict(self.IPMC_OPTS)})
    # nxc is one number for the whole horizon, not one per stage
    with self.assertInException("cannot be cast to OT_INT"):
      self.nxc_solve(p, "auto", [nc]*(N+1))
    # more constant states than states (nx is 2+nc)
    with self.assertInException("is not in [0, nx["):
      self.nxc_solve(p, "auto", 3+nc)
    # negative
    with self.assertInException("is not in [0, nx["):
      self.nxc_solve(p, "auto", -1)

    # ... and a state that simply is not constant.  The state is
    # [p, v, c_0 .. c_{nc-1}], so nxc = nc+1 declares v constant as well,
    # and v's dynamics row is the integrator, not the identity.  Nothing
    # about the declaration is malformed; only the problem disagrees with
    # it, which is exactly the case a range check cannot catch.
    with self.assertInException("declared constant"):
      self.nxc_solve(p, "auto", nc+1)

  @requires_nlpsol("ipmc")
  def test_ipmc_nxc_restoration(self):
    """nxc through the RESTORATION phase, proved rather than hoped for.

    The restoration phase solves a DIFFERENT problem -- minimize the
    infeasibility of yours -- and it computes its steps with the same
    Riccati recursion, so it is the same factorization nxc changes.  It
    was the one part of the algorithm no nxc test reached: every fixture
    here converges from the normal phase and never needs it.

    The trigger is an unreachable terminal position: p_goal = 20 with
    |u| <= 1 over T = 1 from rest, which no trajectory can meet.  ipmc
    stalls in the filter line search, enters restoration, and ends by
    reporting failure -- so this test asserts on a solve that does NOT
    succeed, which is the honest shape of the only fixture that gets
    there.  Both terminal-constraint shapes are used, and they end
    differently: terminal_eq gives return_flag 2 (restoration converged,
    to a point still infeasible for the original problem), the
    inequality form gives return_flag 1 (gave up).  Between them the
    counters say 31 of 40 and 14 of 33 iterations were restoration ones.

    THE PROOF that restoration happened is asserted, not assumed:
    stats()["ipmc"]["restoration_iterations_count"] is nonzero.  Without
    it a fixture that quietly stopped triggering restoration would leave
    this test passing and covering nothing.

    WHAT IS COMPARED.  The sharp assertions are the discrete ones: the
    same total iteration count, the same number of restoration
    iterations, the same return_flag, and bit-identity at nxc = 0.  An
    iterate path this long (33-40 steps, in and out of the restoration
    phase) is an extremely sensitive test of the factorization -- a
    wrong recursion would not stay in step for 40 iterations.  The
    numeric tolerances are the backstop and are much looser than
    NXC_DIGITS, because the final point is not a solution: it is where
    ipmc gave up, with multipliers of magnitude 1e7 to 1e9 that are an
    infeasibility ray rather than a converged dual.  Measured over the
    20 configurations below the worst differences were |df| 3.2e-08,
    |dx| 2.9e-10, |dg| 2.9e-10, |dlam_g| 6.9e-08, |dlam_x| 9.4e-08 in
    checkarray's metric, all of them on the terminal_eq=False pair; the
    terminal_eq=True pair alone would allow four more digits.
    """
    for N, nc, te in [(10, 4, True), (10, 4, False)]:
      p = self.nxc_ocp(N=N, nc=nc, p_goal=20.0, u_max=1.0, terminal_eq=te)
      for sd in ["auto", "manual"]:
        s0, r0 = self.nxc_solve(p, sd, None)
        st0 = s0.stats()
        tag0 = "N=%d nc=%d terminal_eq=%s %s" % (N, nc, te, sd)
        # this fixture is here to reach restoration, and it has to
        self.assertTrue(st0["ipmc"]["restoration_iterations_count"] > 0,
                        tag0+": no restoration iterations, fixture is not "
                             "testing what it says it tests")
        self.assertFalse(st0["success"], tag0)
        for nxc in range(nc+1):
          tag = tag0 + " nxc=%d" % nxc
          s1, r1 = self.nxc_solve(p, sd, nxc)
          st1 = s1.stats()
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
    """The nxc structure has to survive codegen and serialization.

    p.nxc is emitted as a constant and read back by the same runtime
    header, and it is part of the structure fingerprint that decides
    whether a cached solver may be reused -- both of which are silent
    when they are wrong.
    """
    p = self.nxc_ocp(N=5, nc=3, p_goal=1.5)
    solver, r = self.nxc_solve(p, "auto", p["nxc"])
    self.assertTrue(solver.stats()["success"])
    if IPMC_CODEGEN:
      self.check_codegen(solver, p["bounds"], **IPMC_CODEGEN)
    self.check_serialize(solver, p["bounds"])

  @requires_nlpsol("ipmc")
  def test_ipmc_nxc_resolve(self):
    """One solver object, several solves, changing bounds and nxc.

    nxc joins nu/nx/ng/ng_ineq in the fingerprint the runtime compares
    before reusing a cached ipmc solver.  It cannot actually change
    between solves on one object -- it is an nlpsol option -- but the
    fingerprint got one field longer, and an off-by-one there would
    either rebuild every time (slow, invisible) or compare the wrong
    field (wrong, invisible).  Three solves with different bounds on the
    same object pin the behaviour down.
    """
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

  # =====================================================================
  # Soft constraints (nlpsol's 'S' / 's' / 'f_s' slack API) on a
  # STRUCTURED OCP.
  #
  # What is under test here is the interaction between slacks and the
  # Riccati recursion, so every stage-local case is run at all three
  # structure_detection settings.  That is the whole point: under the
  # generic (expand_slacks) de-sugaring the augmented Jacobian is no
  # longer an A B I / C D staircase -- the 2*ns slack columns sit past
  # sum(nx)+sum(nu) and the relaxation rows past sum(ng)+sum(nx) -- so
  # "manual" and "auto" fail at construction outright and only the dense
  # single-stage "none" survives.  Anything that passes below at "manual"
  # or "auto" therefore proves the *native* slack path, and that
  # structure detection still recognises the user-facing OCP through it.
  #
  # Ground truth is recomputed inside every test, never a hard-coded
  # literal, never the same formulation solved twice, and -- since the
  # reference moved to ipopt -- never the same SOLVER twice either:
  #   * a soft optimum is checked against slack_ocp_reference() below, an
  #     independently written augmentation that relaxes EVERY row and so
  #     has a different row order, solved dense BY IPOPT;
  #   * the two degenerate cases are checked against the hard NLP itself;
  #   * the infeasible cases additionally require the hard twin to fail.
  #
  # Why ipopt.  The reference used to be ipmc with
  # structure_detection="none", which made every test ipmc against
  # ipmc: one solver core, two formulations, so a bug in the barrier,
  # the Riccati recursion or the line search that hit both formulations
  # cancelled out.  ipopt is an independent implementation of the same
  # problem class (primal-dual interior point, filter line search), so it
  # cannot cancel that way.  It is NOT built by default in this tree,
  # hence the @requires_nlpsol("ipopt") on every test below that uses the
  # reference; where it is missing those skip rather than silently falling
  # back to a self-comparison.
  #
  # sqpmethod is still deliberately NOT the reference: it returns a wrong
  # objective with success=False on ubs=0 (79.87 instead of 82.38) and on
  # infeasible-hard (78.93, BELOW the true 121.17).  These have two
  # unrelated causes, and neither is the slack layer.  The first is qrqp:
  # it has no anti-cycling rule, so casadi_qrqp_flip and
  # casadi_qrqp_du_index can undo each other on the same index until
  # max_iter, after which it returns a primal-infeasible iterate.  It is
  # a property of the ROW ORDER, not of the formulation -- permuting the
  # rows of one such QP fixes it about 94 times in 100, and the
  # expansion's own order is one of the unlucky ones.  The second is
  # sqpmethod's merit line search stalling; the hand augmentation fails
  # there too.  ipqp segfaults on these problems.
  #
  # Every test carries a non-vacuity guard: the soft optimum has to differ
  # from the hard one by a stated margin (or the hard twin has to fail
  # outright) and a stated number of slacks has to be genuinely active on
  # each side, so a test whose slacks silently stayed at zero -- or whose
  # s_l half stayed at zero while s_u carried the whole relaxation --
  # cannot pass.  slack_activity() counts with a coarse threshold on
  # purpose; see its docstring.
  #
  # Tolerances, re-derived for ipmc-vs-ipopt (the old ones were tuned
  # for ipmc-vs-ipmc and do not transfer).
  #
  # The dominant term is no longer the two formulations' central paths --
  # it is bound_relax_factor.  ipmc, like ipopt at its defaults,
  # widens every bound by bound_relax_factor*max(1,|bound|) = 1e-8 before
  # starting the barrier and never undoes it, so the point it returns
  # sits ~7e-8 OUTSIDE the feasible set (7e-8 = 1e-8 * v_hi, v_hi = 7)
  # and its objective is correspondingly ~2e-6 BELOW the true optimum.
  # Two checks pin this down on the hard N=20 problem: sqpmethod/qrqp
  # (an active-set method, no barrier) and ipopt with
  # bound_relax_factor=0 both return f = 82.377702269 with a constraint
  # violation of 1e-15, while ipmc -- at tol 1e-9 or 1e-13 alike --
  # and ipopt at its default bound_relax_factor both return
  # f = 82.377700303 with a violation of 7.0e-8.
  #
  # The reference therefore runs with bound_relax_factor=0 and solves the
  # problem exactly as written.  Leaving it at the default instead would
  # buy artificially tight agreement on the plain cases (the two solvers
  # would be making the *same* 1e-8 concession) at the price of a much
  # worse artefact on the finite-ubs ones, where it also relaxes the
  # slacks' own 0 <= s <= ubs and moves f by up to 1.1e-5.  Measured both
  # ways; bound_relax_factor=0 wins on every case.
  #
  # Measured worst case over all eighteen problems and all three
  # settings: f 6.8e-6, x 1.8e-6, g 2.2e-7 (all three at ubs_equality_N20,
  # the case with the most active rows), s 8.8e-6 but only for L2, whose
  # should-be-zero slacks sit on the barrier floor sqrt(mu/2Z) of
  # whichever formulation they live in (an L2 penalty has zero gradient at
  # s = 0); every other case has s <= 2.5e-7.  Excluding the finite-ubs
  # cases: f 5.6e-6, x 8.6e-7, g 9.3e-8.
  #
  # Digits are the tightest value leaving >= ~2x headroom on that:
  # digits=4 on f, 5 on x, 6 on g, 4 on s where L2 is in the case list
  # and 5 where it is not.  test_ipmc_slacks_ubs_finite keeps digits=5
  # on g (its 2.2e-7 would leave only 2.3x at digits=6) and gains
  # digits=5 on x, up from the 4 the old reference needed.  The two
  # degenerate tests below compare against the hard NLP rather than
  # against the reference and stay at digits=7.
  # =====================================================================
  IPMC_SLACK_OPTS = {"print_level": 0, "tol": 1e-9, "max_iter": 500}

  # Options for the ipopt reference.  Tight, and bound_relax_factor=0 so
  # that the reference is the optimum of the problem as written; see the
  # tolerance note above.
  IPOPT_REF_OPTS = {"print_time": False,
                    "ipopt": {"print_level": 0, "sb": "yes", "tol": 1e-12,
                              "constr_viol_tol": 1e-12, "dual_inf_tol": 1e-9,
                              "compl_inf_tol": 1e-12, "acceptable_tol": 1e-12,
                              "bound_relax_factor": 0, "max_iter": 3000}}

  def slack_ocp(self, N, soften="band", group_mode="per_row", v_lo=1.0,
                v_hi=7.0, p_goal=5.0, T=1.0, u_max=40.0,
                slope=5.5) -> Dict[str, Any]:
    """A minimum-effort sprint OCP, laid out exactly as ipmc wants.

        p' = v,   v' = u - 0.02 v^2        RK4, N intervals, horizon T

    Decision vector   x = [x_0, u_0, x_1, u_1, ..., x_{N-1}, u_{N-1}, x_N]
    Constraint vector g = [dyn_0, path_0, ..., dyn_{N-1}, path_{N-1}, path_N]

    so nx = [2]*(N+1), nu = [1]*N + [0] and the Jacobian is the staircase.
    The drag term keeps the Hessian non-constant, so the recursion really
    is re-run every iteration.

    Reach p_goal from standstill; the terminal row p_N >= p_goal is always
    hard and is what forces the speed limit to bind.  The speed band on
    stages 1..N is what gets softened -- stage 0's row is left hard,
    because v_0 is pinned by the initial condition and a slack there would
    be vacuous.  That leaves a structurally empty row in S, which the
    slack layer has to tolerate.

    soften      "band"      two-sided path row   v_lo <= v_k <= v_hi
                            (v_lo == v_hi makes it a softened EQUALITY row)
                "upper"     one-sided path row            v_k <= v_hi
                "bound_x"   the speed cap is a simple BOUND -> exercises S_x
                "bound_band_x"  the same, two-sided: v_lo <= v_k <= v_hi as a
                            pair of simple bounds.  The only mode in which the
                            s_l half of an S_x column can be active at all --
                            the sprint starts from standstill, so the early
                            stages sit below the floor while the late ones sit
                            above the ceiling.
                "corridor"  two path rows per stage sharing one slack column
    group_mode  "per_row"   one slack column per softened row  (L1 / L2)
                "single"    one slack column for all of them   (Linf)
                "mixed"     the ODD stages share one column, the even ones
                            keep their own -- so one CROSS-STAGE column (which
                            is lifted into helper states) and N/2 stage-local
                            ones (which stay ipmc slack pairs) in the same
                            problem.  Neither feature's own tests reach that
                            combination.
    """
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
    soft = []   # (row index in the stacked [g; x], group key)

    def group(k):
      if group_mode == "single":
        return "all"
      if group_mode == "mixed":
        return "all" if k % 2 else ("row", k)
      return ("row", k)

    def paths(k):
      rows = []
      if soften == "corridor":
        rows = [(X[k][1], -inf, v_hi), (X[k][0], -inf, slope*k*dt)]
      elif soften == "band":
        rows = [(X[k][1], v_lo if k >= 1 else -inf, v_hi)]
      elif soften == "upper":
        rows = [(X[k][1], -inf, v_hi)]
      n = 0
      for expr, lb, ub in rows:
        if k >= 1:
          soft.append((len(equality), group(k)))
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

    if soften in ["bound_x", "bound_band_x"]:
      for k in range(1, N+1):
        ubx[3*k+1] = v_hi
        if soften == "bound_band_x":
          lbx[3*k+1] = v_lo
        soft.append((ngu+3*k+1, group(k)))

    columns, rows_, cols_ = {}, [], []
    for r, key in soft:
      if key not in columns:
        columns[key] = len(columns)
      rows_.append(r)
      cols_.append(columns[key])
    ns = len(columns)
    S = ca.Sparsity.triplet(ngu+nxu, ns, rows_, cols_)

    return dict(x=x, f=dt*sum(U[k]*U[k] for k in range(N)), g=G,
                nxu=nxu, ngu=ngu, ns=ns, S=S,
                lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg, equality=equality,
                structure=dict(N=N, nx=[2]*(N+1), nu=[1]*N+[0], ng=ng))

  def slack_case(self, penalty, w, ubs=None, par=None,
                 **kwargs) -> Dict[str, Any]:
    """slack_ocp() plus the penalty, so nlp / hard_nlp / bounds are ready.

    'par' is an optional SX parameter that the penalty may depend on; it
    becomes the NLP's 'p' and its value is then supplied per solve (see
    test_ipmc_slacks_parametric_penalty_refused)."""
    p = self.slack_ocp(**kwargs)
    s = ca.SX.sym("s", 2*p["ns"])
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

  def slack_activity(self, s, ns, tol=1e-3):
    """(number of active s_l, number of active s_u).

    The threshold is deliberately coarse. An L2 penalty has zero gradient
    at s = 0, so a slack that ought to be zero stalls on the interior-point
    barrier floor sqrt(mu/2Z) rather than reaching it -- ~5e-6 here. Any
    threshold in between would do; 1e-3 is far above the floor and far
    below every genuinely active slack in these problems (smallest: 0.088).
    """
    sl = [float(s[i]) for i in range(ns)]
    su = [float(s[ns+i]) for i in range(ns)]
    return (sum(1 for v in sl if v > tol), sum(1 for v in su if v > tol))

  def ipmc_slack_solver(self, p, structure_detection):
    opts = {"structure_detection": structure_detection,
            "equality": p["equality"],
            "ipmc": dict(self.IPMC_SLACK_OPTS),
            "S": p["S"]}
    if structure_detection == "manual":
      opts.update(p["structure"])
    return ca.nlpsol("solver", "ipmc", p["nlp"], opts)

  def ipmc_hard_solve(self, p):
    """The hard twin, dense. Returns (result, success)."""
    solver = ca.nlpsol("hard", "ipmc", p["hard_nlp"],
                       {"structure_detection": "none",
                        "ipmc": dict(self.IPMC_SLACK_OPTS)})
    r = solver(**p["bounds"])
    return r, solver.stats()["success"]

  def slack_ocp_reference(self, p, pval=None):
    """Independently written augmentation of the soft problem, solved by ipopt.

    Every row is relaxed here -- the structurally empty rows of S simply
    reproduce the original bound -- so the row order differs from the one
    nlpsol's own de-sugaring produces and agreement with the native path
    cannot be an artefact of a shared layout.

    The solver is IPOPT, not ipmc: a different formulation solved by
    the same core would still cancel a bug in that core's barrier, its
    Riccati recursion or its line search, and this reference exists
    precisely to catch those.  No slack machinery is involved either --
    the slacks are ordinary bounded variables of a dense NLP here."""
    x, g, s = p["x"], p["g"], p["s"]
    nx, ng, ns = p["nxu"], p["ngu"], p["ns"]
    Sd = ca.DM.ones(p["S"])
    S_g, S_x = Sd[:ng, :], Sd[ng:, :]
    s_l, s_u = s[:ns], s[ns:]
    G = ca.vertcat(g+ca.mtimes(S_g, s_l),    # >= lbg
                   g-ca.mtimes(S_g, s_u),    # <= ubg
                   x+ca.mtimes(S_x, s_l),    # >= lbx
                   x-ca.mtimes(S_x, s_u))    # <= ubx
    lbG = ca.vertcat(ca.DM(p["lbg"]), -inf*ca.DM.ones(ng),
                     ca.DM(p["lbx"]), -inf*ca.DM.ones(nx))
    ubG = ca.vertcat(inf*ca.DM.ones(ng), ca.DM(p["ubg"]),
                     inf*ca.DM.ones(nx), ca.DM(p["ubx"]))
    ubs = inf*ca.DM.ones(2*ns) if p["ubs"] is None else ca.DM(p["ubs"])

    rnlp = {"x": ca.vertcat(x, s), "f": p["f"]+p["nlp"]["f_s"], "g": G}
    extra = {}
    if "par" in p:
      rnlp["p"] = p["par"]
      extra["p"] = pval
    solver = ca.nlpsol("reference", "ipopt", rnlp,
                       dict(self.IPOPT_REF_OPTS,
                            ipopt=dict(self.IPOPT_REF_OPTS["ipopt"])))
    r = solver(x0=0,
               lbx=ca.vertcat(-inf*ca.DM.ones(nx), ca.DM.zeros(2*ns)),
               ubx=ca.vertcat(inf*ca.DM.ones(nx), ubs),
               lbg=lbG, ubg=ubG, **extra)
    self.assertTrue(solver.stats()["success"])
    return {"f": r["f"], "x": r["x"][:nx], "s": r["x"][nx:],
            "g": ca.Function("g", [x], [g])(r["x"][:nx])}

  def check_slack_structure(self, solver, p, structure_detection):
    """The OCP structure must survive the slacks untouched.

    With the reference expansion "auto" mis-detects catastrophically --
    it sweeps every slack column into nu[N] and every relaxation row into
    ng[N], collapsing the horizon into the last stage -- so this is what
    tells a working native path from a silently degenerate one."""
    if structure_detection == "none":
      return
    st = solver.stats()
    self.assertEqual(st["N"], p["structure"]["N"])
    for k in ["nx", "nu", "ng"]:
      self.checkarray(ca.DM(p["structure"][k]), ca.DM(st[k]),
                      "structure:"+k, digits=12)

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  def test_ipmc_slacks(self):
    """L1, L2, softened simple bounds (S_x) and a within-stage shared slack.

    Five of the eight stage-local cases; the other three (ubs=0, inactive
    and infeasible-hard) are degenerate and get their own tests.

    l1_bound_band_x_N10 softens a TWO-sided simple bound, which is the only
    case here in which more than one s_l is active at once: the sprint
    starts from standstill, so the first stages sit below the speed floor
    while the last ones sit above the ceiling.  The one-sided cases below
    it (bound_x, corridor) cannot exercise the s_l path at all."""
    L1 = lambda s, w: w*ca.sum1(s)
    L2 = lambda s, w: w*ca.dot(s, s)
    cases = [
      # name, penalty, w, kwargs, (ns, nnz(S)), min df vs hard, min max|s|,
      # (min active s_l, min active s_u)
      ("l1_band_N20",     L1, 0.5, dict(N=20), (20, 20), 1.0, 0.1, (1, 5)),
      ("l2_band_N20",     L2, 2.0, dict(N=20), (20, 20), 1.0, 0.1, (1, 5)),
      ("l1_bound_x_N10",  L1, 0.5, dict(N=10, soften="bound_x", v_hi=6.0),
       (10, 10), 10.0, 1.0, (0, 5)),
      # a two-sided simple bound: BOTH halves of s are active at the optimum
      ("l1_bound_band_x_N10", L1, 0.5,
       dict(N=10, soften="bound_band_x", v_lo=3.0, v_hi=6.0),
       (10, 10), 10.0, 1.0, (2, 5)),
      # two path rows per stage sharing ONE column: the coupling has two
      # entries in the same CD block row range but never leaves the stage
      ("stage_shared_N8", L1, 1.0, dict(N=8, soften="corridor", v_hi=6.0),
       (8, 16), 10.0, 1.0, (0, 4)),
    ]
    for name, penalty, w, kwargs, shape, df_min, s_min, act in cases:
      p = self.slack_case(penalty, w, **kwargs)
      self.assertEqual((p["ns"], p["S"].nnz()), shape)
      ref = self.slack_ocp_reference(p)
      rh, hard_ok = self.ipmc_hard_solve(p)
      self.assertTrue(hard_ok)

      # Non-vacuity: softening has to buy something, and the slacks have
      # to be genuinely active. Without this a test whose slacks all came
      # back zero would pass by reproducing the hard problem.  The counts
      # pin WHICH slacks are active, so a break confined to the s_l path
      # cannot hide behind a large s_u.
      nsl, nsu = self.slack_activity(ref["s"], p["ns"])
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
        # digits: see the tolerance note above.  s stays at 4 here because
        # this is the case list that carries L2 (worst 8.8e-6); the other
        # cases in it are all below 1.4e-7.
        for k, d in [("f", 4), ("x", 5), ("g", 6), ("s", 4)]:
          self.checkarray(ref[k], r[k], name+":"+sd+":"+k, digits=d)
        # feasible for the relaxed bounds, hard rows included
        Sd = ca.DM.ones(p["S"])
        z = ca.vertcat(r["g"], r["x"])
        lb = ca.vertcat(ca.DM(p["lbg"]), ca.DM(p["lbx"]))
        ub = ca.vertcat(ca.DM(p["ubg"]), ca.DM(p["ubx"]))
        self.assertTrue(float(ca.mmax(z-ub-ca.mtimes(Sd, r["s"][p["ns"]:]))) < 1e-6)
        self.assertTrue(float(ca.mmax(lb-ca.mtimes(Sd, r["s"][:p["ns"]])-z)) < 1e-6)
        # The slack description has to survive a round trip: now that the
        # penalty settles at construction, z and Z are VALUES in the stream
        # rather than a Function to re-evaluate, so a solver that lost them
        # would come back with a free slack rather than fail to load.
        self.check_serialize(solver, p["bounds"])

  @requires_nlpsol("ipmc")
  # No @requires_nlpsol("ipopt") here on purpose: this is the one slack test
  # whose ground truth is the hard NLP rather than slack_ocp_reference(), so
  # it still says something in a tree without ipopt.
  def test_ipmc_slacks_degenerate(self):
    """The two cases that must reduce exactly to the hard problem.

    ubs=0 pins every slack; 'inactive' widens the band so far that no
    softened row is active at the optimum.  Both are checked against the
    hard NLP, at all three structure_detection settings -- a native path
    that dropped the slacks on the floor would pass these two and fail
    every other test, and vice versa, which is why they are worth having.
    The non-vacuity guard is the mirror image of the other tests': here
    the slacks must be zero AND the hard problem must be the one that is
    reproduced."""
    L1 = lambda s, w: w*ca.sum1(s)
    cases = [
      ("ubs_zero_N20", dict(N=20), ca.DM.zeros(2*20)),
      ("inactive_N10", dict(N=10, v_lo=-50.0, v_hi=50.0), None),
    ]
    for name, kwargs, ubs in cases:
      p = self.slack_case(L1, 0.5, ubs=ubs, **kwargs)
      self.assertEqual(p["ns"], kwargs["N"])
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
        # The slacks are inert ...
        self.checkarray(ca.DM.zeros(2*p["ns"]), r["s"], name+":"+sd+":s", digits=6)
        # ... so the answer is the hard one. Nothing is being compared
        # across two formulations here, so this is asserted much harder
        # than the tests above: worst observed is 1.4e-13 (ubs=0) and
        # 1.8e-9 (inactive, where the slacks sit on the barrier floor).
        for k in ["f", "x", "g"]:
          self.checkarray(rh[k], r[k], name+":"+sd+":"+k, digits=7)

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  def test_ipmc_slacks_infeasible(self):
    """The headline case: the hard problem is infeasible, the soft one is not.

    Reaching p(T)=5 in T=1 needs a mean speed of 5, but the speed is
    capped at 4, so the hard OCP is infeasible by construction. Softening
    the cap makes it solvable. Non-vacuity here is the strongest of all:
    the hard twin must actually FAIL."""
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
  def test_ipmc_slacks_shared_column(self):
    """Linf: ONE slack column shared by a row of every stage.

    s_u is the largest violation over the whole horizon, so the column is
    a global variable coupling all N stages -- an arrow-head, not a stage
    block.  At structure_detection="none" the whole problem is a single
    stage, so the column is trivially stage-local and ipmc handles it
    exactly.  At "manual"/"auto" it is not: the column has to be LIFTED
    into a helper state held constant along the horizon, which is a
    completely different code path from the stage-local one every other
    test in this file exercises.  Both are run here -- that lifting is
    only reached when a structure is actually detected, so "none" alone
    would leave it untested."""
    L1 = lambda s, w: w*ca.sum1(s)
    cases = [
      ("linf_band_N20", 1.0, dict(N=20, group_mode="single"), 1.0, 0.1, (1, 1)),
      ("linf_infeasible_N12", 5.0,
       dict(N=12, group_mode="single", soften="upper", v_hi=4.0, p_goal=5.0),
       None, 1.0, (0, 1)),
    ]
    for name, w, kwargs, df_min, s_min, act in cases:
      p = self.slack_case(L1, w, **kwargs)
      self.assertEqual(p["ns"], 1)
      ref = self.slack_ocp_reference(p)
      rh, hard_ok = self.ipmc_hard_solve(p)
      nsl, nsu = self.slack_activity(ref["s"], p["ns"])
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
        # The lifted helper state must NOT show up in the reported
        # structure: the user asked for nx=[2]*(N+1) and gets it back.
        self.check_slack_structure(solver, p, sd)
        # ... and the lifting really is what ran, rather than being inferred
        # from the answer coming out right: at "none" nothing is lifted.
        self.assertEqual(solver.stats()["n_lift"], 0 if sd == "none" else 2)
        for k, d in [("f", 4), ("x", 5), ("g", 6), ("s", 5)]:
          self.checkarray(ref[k], r[k], name+":"+sd+":"+k, digits=d)
        # Round trip at every setting, so the LIFTED description -- the helper
        # entities, the rewrite index maps and the per-stage helper row map --
        # is covered as well as the plain slack-pair one at "none".
        self.check_serialize(solver, p["bounds"])

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  def test_ipmc_slacks_lifted_parametric(self):
    """A LIFTED slack column and a PARAMETER in the same problem.

    Every other lifted fixture here has np == 0, which left a hole exactly
    where the two features meet.  On the lifted path the solve loop hands the
    oracle its parameter pointer out of d->nlp, and d->nlp is then the
    REWRITTEN nlp -- helper states appended, rows split and added -- not the
    caller's.  In generated C that second structure is a local built field by
    field rather than the wholesale copy the C++ path makes, so a field it
    forgets is read as whatever was on the stack.  A fixture with no
    parameter cannot notice: nothing dereferences the pointer.

    theta buys effort against slack, so the argmin genuinely moves with it,
    and the two values are asserted to give different answers -- otherwise
    the test would pass just as happily on a parameter that never arrived.
    It goes in f, NOT in f_s: a parametric PENALTY is refused, see
    test_ipmc_slacks_parametric_penalty_refused."""
    L1 = lambda s, w: w*ca.sum1(s)
    q = self.slack_case(L1, 1.0, N=20, group_mode="single")
    self.assertEqual(q["ns"], 1)
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
          # The parameter through generated C, on the lifted path.  This is
          # the combination that was uncovered; without a live theta the
          # generated solver would agree no matter what it was handed.
          if IPMC_CODEGEN:
            self.check_codegen(solver, dict(q["bounds"], p=pval),
                               **IPMC_CODEGEN)
          self.check_serialize(solver, dict(q["bounds"], p=pval))

    # theta has to matter, or none of the above proves it arrived.
    self.assertTrue(float(ca.norm_inf(seen[0]-seen[1])) > 1e-3,
                    "theta did not move the solution")

  @requires_nlpsol("ipmc")
  def test_ipmc_slacks_lifted_distinct_jacobian(self):
    """A LIFTED column on a problem whose Jacobian entries are all DISTINCT.

    The lifted problem's Jacobian is the caller's nonzeros rearranged through
    a construction-time map plus constant +-1 entries.  Every other lifted
    fixture here has near-uniform Jacobian values, so a PERMUTED map -- the
    exact bug a misread Sparsity::triplet mapping produces -- still passes
    them: scattering equal values in the wrong order is invisible.  Here the
    dynamics are nonlinear with distinct coefficients per entry, so a wrong
    map moves the answer and the structure check fires on the first
    Jacobian."""
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
    # one slack column softening the path row of EVERY stage -> cross-stage
    S = ca.DM(ng+nw, 1)
    for k in range(N): S[k*(nx+1)+nx, 0] = 1
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
      if sd != "none": opts["equality"] = eq
      solver = ca.nlpsol("solver", "ipmc", nlp, opts)
      r = solver(**bounds)
      self.assertTrue(solver.stats()["success"], sd)
      self.assertEqual(solver.stats()["n_lift"], 0 if sd == "none" else 2)
      res[sd] = r
      if sd == "auto" and IPMC_CODEGEN:
        self.check_codegen(solver, bounds, **IPMC_CODEGEN)
    # at "none" the whole problem is one stage, the column is stage-local
    # and ipmc's own slack machinery is the (independent) reference
    for k, d in [("f", 5), ("x", 5), ("g", 5), ("s", 5)]:
      self.checkarray(res["none"][k], res["auto"][k], "distinct:"+k, digits=d)

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  def test_ipmc_slacks_mixed_columns(self):
    """A cross-stage column and stage-local ones in the SAME problem.

    This is the first fixture in which ipmc is handed nxc > 0 and ns > 0 at
    once.  Every other Linf problem here shares ONE column, so everything is
    lifted and ipmc's slack machinery is idle; every L1/L2 problem is
    entirely stage-local, so nothing is lifted and nxc stays 0.  Here the
    odd stages share one column -- lifted into two helper states, which are
    declared to ipmc as constant states, hence nxc = 2 -- while each even
    stage keeps its own column and stays an ordinary ipmc slack pair.  The
    two mechanisms then have to coexist inside one factorization: the rank-1
    Schur corrections of the slack pairs are applied to a stage whose state
    block ends in the constant components of the lifted column.

    Non-vacuity is asserted on both halves: the shared column's slack has to
    be active at the optimum AND several of the per-stage ones, or the test
    would be exercising one mechanism with the other switched off in effect.

    Tolerances, measured over the three settings against the ipopt reference:
    f 4.8e-06, x 1.5e-07, g 7.1e-08, s 1.1e-07.  The digits below are the
    usual ones for this block (see the bound_relax_factor comment above) and
    leave 10x on f, 32x on x, 7x on g and 450x on s.
    """
    L1 = lambda s, w: w*ca.sum1(s)
    p = self.slack_case(L1, 5.0, N=12, group_mode="mixed",
                        soften="upper", v_hi=4.0, p_goal=5.0)
    # one shared column (stages 1,3,..,11) + one per even stage
    self.assertEqual(p["ns"], 7)
    ref = self.slack_ocp_reference(p)
    nsl, nsu = self.slack_activity(ref["s"], p["ns"])
    shared_u = float(ref["s"][p["ns"]])       # column 0 is the shared one
    local_u = [float(ref["s"][p["ns"]+j]) for j in range(1, p["ns"])]
    n_local = sum(1 for v in local_u if v > 1e-3)
    print("test_ipmc_slacks_mixed_columns shared", shared_u,
          "active per-stage", n_local, "of", p["ns"]-1)
    self.assertTrue(shared_u > 0.1)           # the lifted column is used
    self.assertTrue(n_local >= 3)             # and so are the slack pairs

    for sd in ["none", "manual", "auto"]:
      solver = self.ipmc_slack_solver(p, sd)
      r = solver(**p["bounds"])
      self.assertTrue(solver.stats()["success"], sd)
      self.check_slack_structure(solver, p, sd)
      # at "none" the single stage makes even the shared column stage-local,
      # so nothing is lifted and this is the pure slack-pair formulation
      self.assertEqual(solver.stats()["n_lift"], 0 if sd == "none" else 2)
      for k, d in [("f", 4), ("x", 5), ("g", 6), ("s", 4)]:
        self.checkarray(ref[k], r[k], "mixed:"+sd+":"+k, digits=d)
      # The only fixture carrying constant states AND ipmc slack pairs at once,
      # so it is the only round trip exercising both halves of the stream.
      self.check_serialize(solver, p["bounds"])

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  def test_ipmc_slacks_lifted_with_nxc(self):
    """A lifted column on top of a user-declared nxc: two trailing blocks.

    The helper states of a lifted column are appended at the END of every
    stage's state and their dynamics are the identity, so they are declared
    to ipmc as constant states.  When the USER has also declared nxc of
    them, the two blocks have to compose into the single contiguous trailing
    block of nxc + n_lift that ipmc's one-number contract describes -- the
    user's constants first, the helpers after them, the same size on every
    stage.  Nothing else in the suite puts both there at once.

    The sharp instrument is the iteration count, not the answer.  Declaring
    nxc does not change the problem at all -- the constant dynamics are
    written out as ordinary gap-closing rows either way -- so a correct
    interface reproduces the iterate path exactly, and a factorization fed
    the wrong block would still converge to the same optimum while taking a
    different number of steps.  What the count cannot see, the independent
    ipopt reference does, and so does the "none" formulation, which reaches
    the same answer through a ipmc slack pair with no lifting at all.

    TOLERANCES, measured over both structure_detection settings.  nxc-on vs
    nxc-off: f 5.3e-15, x 1.1e-15, g 8.9e-16, lam_g 1.7e-13, lam_x 5.3e-15,
    with equal iteration counts (16 and 16) -- NXC_DIGITS leaves 94x on f and
    300x on lam_g.  Against ipopt: f 5.3e-07, x 9.1e-08, g 2.0e-08, s 9.9e-09;
    against the "none" twin: f 9.1e-11, x 2.0e-12, g 2.0e-12, s 1.8e-10.
    """
    q = self.nxc_ocp(N=8, nc=3, p_goal=2.0)
    ns, w = 1, 0.5
    # ONE column softening the speed cap of stages 1..N: cross-stage, so it
    # is lifted, while stage 0's cap row stays hard.
    rows = q["cap_rows"][1:]
    nxu, ngu = q["nlp"]["x"].numel(), q["nlp"]["g"].numel()
    S = ca.Sparsity.triplet(ngu+nxu, ns, rows, [0]*len(rows))
    sv = ca.SX.sym("s", 2*ns)
    p = dict(q, x=q["nlp"]["x"], f=q["nlp"]["f"], g=q["nlp"]["g"],
             nxu=nxu, ngu=ngu, ns=ns, S=S, s=sv, ubs=None,
             nlp={"x": q["nlp"]["x"], "f": q["nlp"]["f"], "g": q["nlp"]["g"],
                  "s": sv, "f_s": w*ca.sum1(sv)})
    ref = self.slack_ocp_reference(p)
    print("test_ipmc_slacks_lifted_with_nxc s", list(numpy.array(ref["s"]).ravel()))
    # non-vacuous: the softened cap really is violated at the optimum
    self.assertTrue(float(ref["s"][ns]) > 0.05)

    def solve(sd, nxc):
      opts = {"structure_detection": sd, "equality": q["equality"],
              "print_time": False, "ipmc": dict(self.IPMC_SLACK_OPTS), "S": S}
      if sd == "manual":
        opts.update(q["structure"])
      if nxc is not None:
        opts["nxc"] = nxc
      solver = ca.nlpsol("solver", "ipmc", p["nlp"], opts)
      return solver, solver(**q["bounds"])

    # "none" is the twin formulation: one dense stage, so the column is
    # stage-local, goes through a ipmc slack pair and nothing is lifted.
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

    # The identity check must still point at the CASADI rows only.  A helper
    # component has no casadi Jacobian row at all, so the check skips them --
    # they are constant by construction, not by declaration.  Declaring one
    # state too many names v, which is NOT constant, and must still be caught:
    # the check has to be reading the right rows with the lifting on.
    with self.assertInException("declared constant"):
      solve("manual", q["nc"]+1)

    # The lifted description through codegen and serialization: n_lift, the
    # lifted partition nxl and the column/side maps are emitted as constants
    # and read back by the same runtime header, and the trailing scalar of
    # the cache fingerprint is now nxc + n_lift -- all silent when wrong.
    sc, _ = solve("auto", q["nxc"])
    if IPMC_CODEGEN:
      self.check_codegen(sc, q["bounds"], **IPMC_CODEGEN)
    self.check_serialize(sc, q["bounds"])

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  def test_ipmc_slacks_ubs_finite(self):
    """A slack that SATURATES at a finite, nonzero upper bound.

    ubs=inf and ubs=0 are both covered elsewhere, and neither reaches this
    active set: here the softened row is still violated at the optimum but
    the slack is not allowed to grow any further, so s sits exactly on ubs
    and the row behaves like a hard row shifted by ubs.  ipmc's own
    barrier carries a `ubsl - sl` / `ubsu - su` term for exactly this, and
    nothing was exercising it.

    Non-vacuity is the saturation count itself, asserted per case: a run
    whose slacks all stayed strictly inside their bounds would prove
    nothing over the ubs=inf tests.  The bounds are deliberately DIFFERENT
    on the two halves, so an implementation that fed ubsl where ubsu
    belongs would land on a different (or infeasible) problem.

    Tolerances. These three cases carry the worst discrepancy in the file
    -- the optimum is a vertex at which up to eleven slacks are pinned
    against ubs on top of the usual active rows, so ipmc's 1e-8 bound
    relaxation (see the block comment above) is applied to more active
    rows here than anywhere else.  Measured worst case over the three
    cases and all three settings, against the ipopt reference: f 6.8e-6,
    x 1.8e-6, g 2.2e-7, s 2.5e-7 (all four at ubs_equality_N20).  x gains
    a digit over the ipmc-vs-ipmc era (5, not 4) and g gives one up
    (5, not 6): at digits=6 its 2.2e-7 would leave only 2.3x headroom.
    """
    L1 = lambda s, w: w*ca.sum1(s)
    cases = [
      # name, w, kwargs, (ubsl, ubsu), min (saturated s_l, saturated s_u),
      # hard twin feasible?
      # a two-sided SIMPLE BOUND (S_x): both halves saturate
      ("ubs_bound_band_x_N10", 0.5,
       dict(N=10, soften="bound_band_x", v_lo=3.0, v_hi=6.0),
       (1.0, 0.6), (1, 2), True),
      # softened equality rows, capped hard: eleven of the twenty upper
      # slacks pile up on ubs
      ("ubs_equality_N20", 0.5, dict(N=20, v_lo=3.0, v_hi=3.0),
       (1.8, 3.0), (1, 11), False),
      # the same on a LIFTED (cross-stage, Linf) column: ubs is carried by
      # the helper state's own simple bound, not by a ipmc slack pair
      ("ubs_linf_N20", 1.0, dict(N=20, group_mode="single"),
       (0.1, 0.1), (1, 1), True),
    ]
    for name, w, kwargs, cap, sat_min, hard_feasible in cases:
      ns = 1 if kwargs.get("group_mode") == "single" else kwargs["N"]
      ubs = ca.vertcat(cap[0]*ca.DM.ones(ns), cap[1]*ca.DM.ones(ns))
      p = self.slack_case(L1, w, ubs=ubs, **kwargs)
      self.assertEqual(p["ns"], ns)
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
      # ... and they are AT the bound, not past it. 1e-6 rather than 0:
      # ipmc relaxes every bound by a small factor before starting the
      # barrier, so a saturated slack lands ~3e-8 above ubs, in both
      # formulations. Same tolerance as the relaxed-feasibility checks in
      # test_ipmc_slacks.
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
        for k, d in [("f", 4), ("x", 5), ("g", 5), ("s", 5)]:
          self.checkarray(ref[k], r[k], name+":"+sd+":"+k, digits=d)
        # ubs is a runtime INPUT, not part of the description, so the round
        # trip is fed it explicitly -- a deserialised solver that dropped the
        # cap would saturate elsewhere and the digits=16 compare says so.
        self.check_serialize(solver, dict(p["bounds"], ubs=ubs))

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  def test_ipmc_slacks_equality_row(self):
    """A softened EQUALITY row (lbg == ubg).

    This is a silent-failure trap rather than an exotic feature.  The row
    is tagged as an equality in the 'equality' option -- structure
    detection needs that -- but once relaxed it is genuinely two-sided,
    lbg - s_l <= g <= ubg + s_u, and must be classified as an INEQUALITY
    for the solver: an equality row has no slack side to attach s to, and
    a path that trusted the 'equality' tag would either drop the
    relaxation or pin the row.  Both failure modes are silent -- the
    solver still converges, to the hard problem's answer or to nothing at
    all -- so only comparing against an independent relaxation catches it.

    v is pinned to 3 at every stage but the terminal row needs a mean
    speed of 5, so the hard twin is infeasible by construction (that is
    the non-vacuity guard) and BOTH halves of s are active at the
    optimum: the first stages cannot accelerate to 3 in one step, the
    later ones have to exceed it."""
    p = self.slack_case(lambda s, w: w*ca.sum1(s), 0.5, N=20,
                        v_lo=3.0, v_hi=3.0)
    # every softened row really is a two-sided equality row, and is
    # declared as one
    for r in set(p["S"].row()):
      self.assertTrue(r < p["ngu"])
      self.assertEqual(p["lbg"][r], p["ubg"][r])
      self.assertTrue(p["equality"][r])

    ref = self.slack_ocp_reference(p)
    rh, hard_ok = self.ipmc_hard_solve(p)
    self.assertFalse(hard_ok)
    nsl, nsu = self.slack_activity(ref["s"], p["ns"])
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
    """Different weights on the lower and the upper half of s.

    Every other test in this file weights s_l and s_u identically
    (w*sum1(s) or w*dot(s,s)), which makes a zl/zu -- or an S_l/S_u --
    swap invisible.  Here the two halves cost 0.3 and 3.0, on a two-sided
    simple bound where both halves are genuinely active, and the case is
    run BOTH ways round.  Agreement with the reference already pins which
    weight went where; the explicit comparison at the end additionally
    rules out an implementation that averaged them or read one twice."""
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
      ns = p["ns"]
      nsl, nsu = self.slack_activity(ref["s"], ns)
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

    # Swapping the two weights is a different problem, and in the
    # direction the penalty predicts: the run in which the UPPER slack is
    # the cheap one (0.3) uses far more of it.  Measured 6.82 vs 2.21.
    self.assertTrue(sums[(3.0, 0.3)][1] > sums[(0.3, 3.0)][1] + 1.0)
    self.assertTrue(abs(sums[(3.0, 0.3)][2] - sums[(0.3, 3.0)][2]) > 1.0)

  @requires_nlpsol("ipmc")
  @requires_nlpsol("ipopt")
  def test_ipmc_slacks_parametric_penalty_refused(self):
    """A penalty whose weight is the NLP parameter p.

    z and Z are evaluated ONCE, by IpmcInterface::init, and handed to the
    runtime (and to the symbolic rewrite that lifts a cross-stage column)
    as numbers; nothing at solve time evaluates f_s any more.  A penalty
    that depends on p therefore cannot be honoured -- every solve would
    silently use the weights p happened to hold at construction -- so it
    is refused where the fact is knowable, at construction, rather than
    solved wrongly.  That refusal is what this test pins down: an
    implementation that quietly baked in the weights and kept solving
    would pass the rest of this file.

    The capability itself is not lost, it moves: expand_slacks=True
    carries f_s into the NLP, where p reaches it like any other parameter.
    So the second half solves the very same problem through that route,
    for three solves w = 0.3, 3.0 and 0.3 again.  The first two must
    differ in the direction the penalty predicts (a ten-fold dearer slack
    buys less relaxation) and each must match its own reference; the third
    must reproduce the first exactly, which is what a stale cache cannot
    do."""
    par = ca.SX.sym("w_slack")
    p = self.slack_case(lambda s, w: par*ca.sum1(s), None, par=par,
                        N=20, soften="bound_band_x", v_lo=3.0, v_hi=6.0)
    weights = [0.3, 3.0]
    ref = {pv: self.slack_ocp_reference(p, pv) for pv in weights}
    tot = {pv: float(ca.sum1(ref[pv]["s"])) for pv in weights}
    print("test_ipmc_slacks_parametric_penalty_refused sum|s|", tot)
    # Non-vacuity: the two parameter values really do give different
    # answers, and the cheap one relaxes more. Measured 12.47 vs 3.85.
    self.assertTrue(tot[0.3] > tot[3.0] + 1.0)
    self.assertTrue(float(ca.norm_inf(ref[0.3]["x"]-ref[3.0]["x"])) > 0.1)
    for pv in weights:
      nsl, nsu = self.slack_activity(ref[pv]["s"], p["ns"])
      self.assertTrue(nsl >= 1 and nsu >= 9)

    # Refused at construction, on every structure detection, and the
    # message has to say why and what to do instead.
    for sd in ["none", "manual", "auto"]:
      print("test_ipmc_slacks_parametric_penalty_refused", sd)
      with self.assertRaises(Exception) as cm:
        self.ipmc_slack_solver(p, sd)
      msg = str(cm.exception)
      self.assertTrue("depends on the parameter p" in msg, msg)
      self.assertTrue("expand_slacks" in msg, msg)

    # ... and still solvable through the reference expansion, which does
    # track p.
    opts = {"structure_detection": "none", "equality": p["equality"],
            "ipmc": dict(self.IPMC_SLACK_OPTS), "S": p["S"],
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
    # Same p, same solver, same answer -- bit for bit, not to a
    # tolerance: nothing about the second solve may leak into the third.
    for k in ["f", "x", "g", "s"]:
      self.checkarray(out[0][k], out[2][k],
                      "parametric_repeat:expand:"+k, digits=12)

if __name__ == '__main__':
    unittest.main()
