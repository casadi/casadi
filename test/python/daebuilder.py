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
from casadi import *

import casadi as ca

import subprocess
import os
import unittest
from helpers import *

class Daebuildertests(casadiTestCase):

  def test_reference_fmus(self):
    if "ghc-filesystem" not in CasadiMeta.feature_list(): return
    for name in ["VanDerPol2","VanDerPol3"]:
        fmu_file = "../data/" + name + ".fmu"
        if not os.path.exists(fmu_file):
            print("Skipping test_reference_fmus, resource not available")
            return
        dae = DaeBuilder("car",fmu_file)
        dae.disp(True)
        x0 = SX.sym('x0')
        x1 = SX.sym('x1')
        mu = SX.sym("mu")
        
        x = vertcat(x0,x1)
        c = mu
        f = dae.create('f',['x'],['ode'])
        f_ref = Function('f',[x],[vertcat(x1,1 * ((1 - x0 * x0) * x1) - x0)])
        test_point = [vertcat(1.1,1.3)]
        self.checkfunction(f,f_ref,inputs=[vertcat(1.1,1.3)],digits=7)
        if not name.endswith("3"):
            self.check_serialize(f,inputs=test_point)
  
  def test_cstr(self):
    fmu_file = "../data/cstr.fmu"
    if not os.path.exists(fmu_file):
        print("Skipping test_fmu_zip, resource not available")
        return
    unzipped_name = "cstr"
    unzipped_path = os.path.join(os.getcwd(), unzipped_name)
    import shutil
    if os.path.isdir(unzipped_path): shutil.rmtree(unzipped_path)
    import zipfile
    with zipfile.ZipFile(fmu_file, 'r') as zip_ref:
        zip_ref.extractall(unzipped_name)
    dae = DaeBuilder("cstr",unzipped_name)
    f = dae.create('f',['x','u'],['ode'])
    with self.assertInException("No stats available"):
        f.stats()
    dae = None
    
    C_A = SX.sym('C_A')
    C_B = SX.sym('C_B')
    V = 1
    k = 0.5
    
    q_in = SX.sym('q_in')
    C_A_in = SX.sym('C_A_in')
    C_B_in = SX.sym('C_B_in')

    x = vertcat(C_A,C_B)
    u = vertcat(C_A_in,C_B_in,q_in)    
    ode = vertcat((q_in*(C_A_in - C_A) - V*k*C_A*C_B)/V,(q_in*(C_B_in - C_B) + V*k*C_A*C_B)/V)


    f_ref = Function('f',[x,u],[ode],['x','u'],['ode'])
    
    test_point = [vertcat(1.1,1.3),vertcat(1.7,1.11,1.13)]
    self.checkfunction(f,f_ref,inputs=test_point,digits=4,hessian=False,evals=1)
    print(f.stats())
          
  def test_rumoca(self):
    if "rumoca" not in CasadiMeta.feature_list(): return
    rumoca = os.path.join(GlobalOptions.getCasadiPath(),'rumoca')
    rumoca_exe = os.path.join(GlobalOptions.getCasadiPath(),'rumoca.exe')
    if os.path.exists(rumoca) or os.path.exists(rumoca_exe):
        p = subprocess.run([rumoca,"-t","../assets/casadi_daebuilder.jinja","-m", "../assets/hello_world.mo"]) 

  def test_fmu_zip(self):
    fmu_file = "../data/VanDerPol2.fmu"
    if not os.path.exists(fmu_file):
        print("Skipping test_fmu_zip, resource not available")
        return
    unzipped_name = "VanDerPol2"
    unzipped_path = os.path.join(os.getcwd(), unzipped_name)
    for serialize_mode in ["link","embed"]:
        for use_zip in [True]:#,False]:
            if use_zip:
                if "ghc-filesystem" in CasadiMeta.feature_list():
                    dae = DaeBuilder("cstr",fmu_file,{"resource_serialize_mode": serialize_mode})
                else:
                    with self.assertInException("passing fmu files to DaeBuilder is unsupported"):
                        dae = DaeBuilder("cstr",fmu_file)
                    continue
            else:
                import shutil
                if os.path.isdir(unzipped_path): shutil.rmtree(unzipped_path)
                import zipfile
                with zipfile.ZipFile(fmu_file, 'r') as zip_ref:
                    zip_ref.extractall(unzipped_name)
                print('Unzipped %s into %s' % (fmu_file, unzipped_path))
                dae = DaeBuilder("cstr",unzipped_path)
            dae.disp(True)
            f = dae.create('f',['x'],['ode'])
            dae = None
            
            x0 = SX.sym('x0')
            x1 = SX.sym('x1')
            mu = SX.sym("mu")
            
            x = vertcat(x0,x1)
            c = mu
            f_ref = Function('f',[x],[vertcat(x1,1 * ((1 - x0 * x0) * x1) - x0)])
            test_point = [vertcat(1.1,1.3)]

            self.checkfunction_light(f,f_ref,inputs=test_point,digits=7)
            
            f.save('f.casadi')
            
            f = None
            
            f = Function.load("f.casadi")
            self.checkfunction_light(f,f_ref,inputs=test_point,digits=7)
            
            # Type decay, so test twice
            f.save('f.casadi')
            self.checkfunction_light(f,f_ref,inputs=test_point,digits=7)
            f = None
            f = Function.load("f.casadi")
            self.checkfunction(f,f_ref,inputs=test_point,hessian=False,digits=7,evals=1)
  
  @requires_rootfinder("kinsol")
  @requires_rootfinder("newton")
  @requires_nlpsol("ipopt")
  @requires_integrator("cvodes")
  def test_fmu_demo(self):
        # Use FMU to create a CasADi/DaeBuilder instance

        from time import time

        t1 = time()
        fmu_file = '../data/vdp.fmu'
        if not os.path.exists(fmu_file):
            print("Skipping test_fmu_demo, resource not available")
            return
        dae = ca.DaeBuilder('vdp', fmu_file)
        dae.disp(True)

        print('FMU ({}) provides directional derivatives: {}'.format(fmu_file, dae.provides_directional_derivative()))

        # Get state vector, initial conditions, bounds
        x = dae.x()
        lbx = dae.min(x)
        ubx = dae.max(x)
        x0 = dae.start(x)
        print('x: ', x)
        print('lbx: ', lbx)
        print('ubx: ', ubx)
        print('x0: ', x0)

        # Get free control, bounds
        u = dae.u()
        lbu = dae.min(u)
        ubu = dae.max(u)
        print('u: ', u)
        print('lbu: ', lbu)
        print('ubu: ', ubu)

        # Evaluate ODE right-hand-side
        f = dae.create('f', ['x', 'u'], ['ode'])
        print(f)
        # Evaluate the function numerically
        u_test = 0.4
        xdot_test = f(x0, u_test)
        print('xdot_test: ', xdot_test)



        # Alternative ODE right-hand-side function with u fixed
        dae.set('u', 0.4)
        f_no_u = dae.create('f_no_u', ['x'], ['ode'])
        print(f_no_u)
        # Evaluate the function numerically
        xdot_test_no_u = f_no_u(x0)
        print('xdot_test_no_u: ', xdot_test_no_u)

        # Evaluate ODE right-hand-side with auxilliary field
        f_with_aux = dae.create('f_with_aux', ['x', 'u'], ['ode'],
                       dict(aux = ['x1', 'x2', 'u']))
        print(f_with_aux)
        # After evaluation, the aux variables are in the stats field:
        u_test = 0.4
        xdot_test_with_aux = f_with_aux(x0, u_test)
        print('xdot_test_with_aux: ', xdot_test_with_aux)
        print(f_with_aux.stats())

        # A function that calculates the ODE and two Jacobian blocks:
        J = dae.create('J', ['x', 'u'], ['ode', 'jac_ode_x', 'jac_ode_u'],
                      dict(verbose = True, parallelization = 'openmp'))
        print(J)

        # Evaluate Jacobian
        xdot_test, A_test, B_test = J(x0, u_test)
        print('xdot_test: ', xdot_test)
        print('A_test: ', A_test)
        print('B_test: ', B_test)

        # A function that calculates the ODE and two Jacobian blocks:
        J_fd = dae.create('J_fd', ['x', 'u'], ['ode', 'jac_ode_x', 'jac_ode_u'],
                      dict(verbose = True,
                           parallelization = 'openmp',
                           enable_ad = False))
        print(J_fd)

        # Evaluate Jacobian
        xdot_fd, A_fd, B_fd = J_fd(x0, u_test)
        print('xdot_fd: ', xdot_fd)
        print('A_fd: ', A_fd)
        print('B_fd: ', B_fd)


        # Create a function that calculate the Hessian of a linear
        # combination of "ode" with respect to x and u: 
        H = dae.create('H', ['x', 'u', 'adj_ode'],
                       ['jac_adj_x_x', 'jac_adj_x_u', 'jac_adj_u_u'],
                       dict(verbose = True, parallelization = 'openmp'))
        print(H)
        # Evaluate Hessian
        # Evaluate Jacobian
        H_xx, H_xu, H_uu = H(x0, u_test, 1)
        print('H_xx: ', H_xx)
        print('H_xu: ', H_xu)
        print('H_uu: ', H_uu)



        import numpy as np
        # Create an integrator for simulating the ODE over 10 s using CVODES
        T = 10.
        tgrid = np.linspace(0., T, 100)
        dae.set_min(x[0], -np.inf) # relax lower bound on x1
        daefun = dae.create('daefun')
        sim = ca.integrator('sim', 'cvodes', daefun, 0, tgrid)




        # Call integrator instance
        r = sim(x0 = x0, u = u_test)
        x_traj_sim = r['xf'].full()



        # Create a new Simulator with support for forward sensitivities
        fwd_sim = ca.integrator('fwd_sim', 'cvodes', daefun, 0, tgrid, dict(nfwd = 1, verbose = False, print_stats = True, second_order_correction = False))
        print(fwd_sim)
        # Let's calculate sensitivity w.r.t. x1, i.e. use a forward seed [1, 0] for x0:
        seed_x0 = ca.DM([1, 0])
        fwd_test = fwd_sim(x0 = ca.horzcat(x0, seed_x0),
                          u = ca.horzcat(u_test, 0))
        # The directional derivatives can be found in column 1, 3, 5, etc. of `xf` for the different grid points.
        # (column 0, 2, 4, etc corresponds to nondifferentiated xf, for each grid point)
        print('dxf/dx0 * seed_x0 [AD, forward] = ', fwd_test['xf'][:, 1::2])
        # We can compare this result with a finite difference perturbation
        pert = 1e-3
        pert_sim = sim(x0=x0 + pert * seed_x0, u=u_test)
        x_pert_sim = pert_sim['xf']
        print('dxf/dx0 * seed_x0 [FD] = ', (x_pert_sim - x_traj_sim) / pert)



        # Create a new Simulator with support for adjoint sensitivities
        adj_sim = ca.integrator('adj_sim', 'cvodes', daefun, 0, tgrid, dict(nadj = 1, verbose = False, print_stats = True, second_order_correction = True, print_time = True))
        print(adj_sim)
        # Let's calculate sensitivity of xf[0] w.r.t. u at all the grid points, i.e. use an adjoint seed [1, 0] for xf:
        seed_xf = ca.DM.zeros(2, 100)
        seed_xf[:, -1] =  ca.DM([1, 0])  # We seed the state at the last grid point only
        adj_test = adj_sim(x0 = x0, adj_xf = seed_xf, u = u_test)
        print('trans(dxf/du) * seed_xf (AD, adjoint) = ', adj_test['adj_u'])

        fwd_adj_sim = ca.integrator('fwd_adj_sim', 'cvodes', daefun, 0, tgrid, dict(nfwd = 1, nadj = 1, verbose = False, print_stats = True, second_order_correction = False))
        print(fwd_adj_sim)
        # Let's calculate sensitivity of d(xf[0])/dx1 w.r.t. u at all the grid points, i.e. combining the seeds for forward and adjoint sensitivitiy analysis above
        fwd_adj_test = fwd_adj_sim(x0 = ca.horzcat(x0, seed_x0), adj_xf = ca.horzcat(seed_xf, ca.DM.zeros(seed_xf.shape)), u = ca.horzcat(u_test, 0))
        print('d(trans(dxf/du) * seed_xf)/dx0 * seed_x0 (AD, forward-over-adjoint) = ', fwd_adj_test['adj_u'][:,1::2])

        # Number of integrator steps
        N = 50
        # Size of the finite elements
        h = T/N

        # Let us use a 4th order collocation discretization using Legendre roots, cf. 
        # Nonlinear Programming: Concepts, Algorithms, and Applications to Chemical Processes
        # by Lorenz Biegler (2010)
        d = 3
        # The roots can be queried from CasADi or looked up in the above textbook
        tau_root = np.append(0, ca.collocation_points(d, 'legendre'))
        # Coefficients of the collocation equation
        C = np.zeros((d+1,d+1))
        # Coefficients of the continuity equation
        D = np.zeros(d+1)
        # Coefficients of the quadrature function
        B = np.zeros(d+1)
        # Construct polynomial basis
        for j in range(d+1):
            # Construct Lagrange polynomials
            p = np.poly1d([1])
            for r in range(d+1):
                if r != j:
                    p *= np.poly1d([1, -tau_root[r]]) / (tau_root[j]-tau_root[r])
            # Evaluate at the final time to get the coefficients of the continuity equation
            D[j] = p(1.0)
            # Evaluate the time derivative to get the coefficients of the continuity equation
            pder = np.polyder(p)
            for r in range(d+1): C[j,r] = pder(tau_root[r])
            # Evaluate the integral of the polynomial to get the coefficients of the quadrature function
            pint = np.polyint(p)
            B[j] = pint(1.0)
        print('tau_root: ', tau_root)
        print('C: ', C)
        print('D: ', D)
        print('B: ', B)

        # Symbolic expression for the controls (piecewise constant)
        U = ca.MX.sym('U',len(u))
        # Symbolic expressions for the states at each collocation point
        X = [ca.MX.sym('X_' + str(j), len(x)) for j in range(d+1)]
        # Define the collocation equations
        g = []
        for j in range(1,d+1):
          # Expression for the state derivative at the collocation point
          xdot_j = 0
          for r in range (d+1): xdot_j += C[r,j]*X[r]
          # Append collocation equations
          g.append(h*f(X[j],U) - xdot_j)
        # Concatenate constraints
        g = ca.vertcat(*g)
        # Form a root-finding function, implicitly defining X[1:] as a function of U, X[0]
        X_unknown = ca.vertcat(*X[1:])
        rfp = ca.Function('rfp', [X_unknown, X[0], U], [g], ['V', 'X0', 'U'], ['res'])

        # We can solve this system of equations using Sundials/KINSOL
        ifcn = ca.rootfinder('ifcn', 'kinsol', rfp, dict(print_level = 1))
        # Take a single step of the integrator
        v0 = ca.repmat(x0, d, 1)
        v0 = ifcn(v0, x0, u_test)
        # State at each collocation point
        x_all = ca.reshape(v0, len(x), d)
        # Prepend initial state
        x_all = ca.horzcat(x0, x_all)
        # Also calculate the state at the end
        xf = ca.mtimes(x_all, D)
        # Print solution
        print('x_all: ', x_all)
        print('xf: ', xf)

        # We can also use CasADi's native Newton solver
        ifcn = ca.rootfinder('ifcn', 'newton', rfp, dict(print_iteration = True))
        # Take a single step of the integrator
        v0 = ca.repmat(x0, d, 1)
        v0 = ifcn(v0, x0, u_test)
        # State at each collocation point
        x_all = ca.reshape(v0, len(x), d)
        # Prepend initial state
        x_all = ca.horzcat(x0, x_all)
        # Also calculate the state at the end
        xf = ca.mtimes(x_all, D)
        # Print solution
        print('x_all: ', x_all)
        print('xf: ', xf)

        # We can solve this system of equations using Sundials/KINSOL
        print(rfp)
        ifcn = ca.rootfinder('ifcn', 'kinsol', rfp, dict(print_level = 1))
        print(ifcn)
        # Differentiate Newton solver
        jac_ifcn = ifcn.jacobian()
        sol = jac_ifcn(V0 = v0, X0 = x0, U = u_test)
        # Jacobian of state at collocation points w.r.t. U
        print('jac_V_U = ', sol['jac_V_U'])
        # Jacobian of state at collocation points w.r.t. X0
        print('jac_V_X0 = ', sol['jac_V_X0'])

        # Compare with finite differences
        pert = 1e-3
        # Perturn u
        v0_pert = ifcn(ca.repmat(x0, d, 1), x0, u_test + pert)
        jac_V_U_fd = (v0_pert - v0) / pert
        print('jac_V_U (FD) = ', jac_V_U_fd)
        # Perturb x1
        v0_pert = ifcn(ca.repmat(x0, d, 1), x0 + ca.DM([pert, 0]), u_test)
        jac_V_X0_fd1 = (v0_pert - v0) / pert
        # Perturb x2
        v0_pert = ifcn(ca.repmat(x0, d, 1), x0 + ca.DM([0, pert]), u_test)
        jac_V_X0_fd2 = (v0_pert - v0) / pert
        print('jac_V_X0 (FD) = ', ca.horzcat(jac_V_X0_fd1, jac_V_X0_fd2))

        # Let's create a CasADi function for simulating the whole trajectory using Kinsol
        ifcn = ca.rootfinder('ifcn', 'kinsol', rfp)
        x0_in = ca.MX.sym('x0', len(x))
        u_all = []
        xk = x0_in
        x_all = [x0_in]
        for k in range(N):
            # Symbolic expression for the control for the interval
            uk = ca.MX.sym('u' + str(k), len(u))
            u_all.append(uk)
            # Solve rootfinding problem, using previous step as guess
            vk = ifcn(ca.repmat(xk,d,1), xk, uk)
            # Reshape and prepend initial state
            vk = ca.reshape(vk, len(x), d)
            vk = ca.horzcat(xk, vk)
            # Get the state at the end
            xk = ca.mtimes(vk, D)
            # Save trajectory
            x_all.append(xk)
        # Embed in a CasADi Function
        irksim = ca.Function('irksim', [x0_in, ca.vcat(u_all)], [ca.hcat(x_all)],
                                 ['x0', 'u'], ['x'])
        # Evaluate the function
        x_traj = irksim(x0, u_test).full()

        # A nonlinear equation solver is a differentiable object in CasADi.
        # We can differentiate it analytically to get the Jacobian-times-vector product:
        fwd_irksim = irksim.forward(1)
        # Let's calculate sensitivity w.r.t. x1, i.e. use a forward seed [1, 0] for x0:
        seed_x0 = ca.DM([1, 0])
        fwd_test = fwd_irksim(x0 = x0, u = u_test, fwd_x0 = seed_x0)
        print('fwd_x (AD, forward) = ', fwd_test['fwd_x'])
        # We can compare this result with a finite difference perturbation
        pert = 1e-3
        x_pert = irksim(x0 + pert * seed_x0, u_test).full()
        print('fwd_x (FD) = ', (x_pert - x_traj) / pert)

        # Now let's to a reverse mode AD, i.e. a transposed-Jacobian-times-vector product
        adj_irksim = irksim.reverse(1)
        # We seed the last entry in x2
        seed_x = ca.DM.zeros(2, N + 1)
        seed_x[1, N] = 1
        adj_test = adj_irksim(x0 = x0, u = u_test, adj_x = seed_x)
        print('adj_x0 (AD, reverse) = ', adj_test['adj_x0'])
        print('adj_u (AD, reverse) = ', adj_test['adj_u'])

        # We can check that the result is consistent with forward mode AD:
        print('Forward mode: derivative of x[1, N] w.r.t. x0[0] = ', fwd_test['fwd_x'][1, N])
        print('Reverse mode: derivative of x[1, N] w.r.t. x0[0] = ', adj_test['adj_x0'][0])

        # We can also calculate the full Jacobian of x w.r.t. x0 and u:
        jac_irksim = irksim.jacobian()
        jac_test = jac_irksim(x0 = x0, u = u_test)
        print('jac_x_x0 = ', jac_test['jac_x_x0'])
        print('jac_x_u = ', jac_test['jac_x_u'])

        # We will construct an NLP incrementally, starting with no decision variables, 
        # and constraints, and zero objective
        w = []
        g = []
        J = 0
        # Loop over all control intervals, constructing the NLP and x, u trajectories
        U_all = []
        X_all = []
        xk = ca.MX(x0)
        for k in range(N):
            # Symbolic expression for the control for the interval
            uk = ca.MX.sym('u' + str(k), len(u))
            U_all.append(uk)
            w.append(uk)
            # States at each collocation point for the interval
            Xc = [xk]
            for j in range(1, d + 1):
                xkj = ca.MX.sym('x' + str(k) + '_' + str(j), len(x))
                Xc.append(xkj)
                X_all.append(xkj)
                w.append(xkj)
            # Collect the collocation equations
            for j in range(1,d+1):
              # Expression for the state derivative at the collocation point
              xdot_j = 0
              for r in range (d+1): xdot_j += C[r,j]*Xc[r]
              # Append collocation equations
              g.append(h*f(Xc[j],uk) - xdot_j)
            # State at the end of the interval
            xk_end = 0
            for j in range(0,d+1): xk_end = xk_end + D[j] * Xc[j]
            # New variable for the state at the end of the interval
            xk = ca.MX.sym('x' + str(k + 1), len(x))
            X_all.append(xk)
            w.append(xk)
            # Enforce continuity
            g.append(xk_end - xk)
            # Add contributions to the objective
            J = J + ca.sumsqr(uk) + ca.sumsqr(xk)
        # Concatenate vectors
        w = ca.vcat(w)
        g = ca.vcat(g)
        X_all = ca.hcat(X_all)
        U_all = ca.hcat(U_all)
        # Create mappings from (u, x) -> (w) and back
        to_w = ca.Function('to_w', [X_all, U_all], [w], ['x', 'u'], ['w'])
        from_w = ca.Function('from_w', [w], [X_all, U_all], ['w'], ['x', 'u'])
        # Form the NLP
        nlp = dict(x = w, f = J, g = g)


        # Let's form the Jacobian of g w.r.t. w
        t0 = time()
        jac_g_w = ca.jacobian(g, w)
        print('Symbolic expression for Jacobian of g w.r.t w formed in %g s' % (time() - t0))
        # Create a function for evaluation
        t0 = time()
        jfcn = ca.Function('jfcn', [w], [jac_g_w], ['w'], ['jac_g_w'])
        print('Jacobian function formed in %g s' % (time() - t0))
        # Evaluate function
        w_test = to_w(u = u_test, x = x0)['w']
        t0 = time()
        J_test = jfcn(w_test)
        print('Jacobian evaluated in %g s' % (time() - t0))

        # Create an IPOPT instance, using L-BFGS
        opts = dict(ipopt = dict(hessian_approximation = 'exact', linear_solver = 'mumps'))
        solver = ca.nlpsol('solver', 'ipopt', nlp, opts)



        # Generate initial guess for w
        w0 = to_w(u = 0, x = x0)['w']
        # Lower bound for w
        lbw = to_w(u = lbu, x = lbx)['w']
        # Upper bound for w
        ubw = to_w(u = ubu, x = ubx)['w']
        # Solve NLP
        sol = solver(x0 = w0, lbx = lbw, ubx = ubw, lbg = 0, ubg = 0)


if __name__ == '__main__':
    unittest.main()
