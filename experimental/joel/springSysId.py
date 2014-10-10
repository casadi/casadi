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
import matplotlib.pyplot as plt
import numpy as N
import casadi as C

raise Exception("This example does not work, INTEGRATOR_T0 and INTEGRATOR_TF must be eliminated")


# Construct an explicit rk4 integrator
def create_integrator_rk4():
    u = C.SX("u") # control for one segment
    k = C.SX("k") # spring constant
    b = C.SX("b") # viscous damping coefficient
  
    # initial state
    s0 = C.SX("s0") # initial position
    v0 = C.SX("v0") # initial speed
  
    t0 = C.SX("t0") # initial time
    tf = C.SX("tf") # final time
  
    m = 1.0   # mass

    # apply fixed step rk4
    def eoms(x_, scalar, kn, t_):
        s_ = x_[0] + scalar*kn[0]
        v_ = x_[1] + scalar*kn[1]

        return [v_, (u - k*s_ - b*v_)/m]
    
    dt = tf-t0 # time step
    
    x0 = [s0, v0]

    k1 = eoms( x0,      0, [0,0], t0          );
    k2 = eoms( x0, 0.5*dt,    k1, t0 + 0.5*dt );
    k3 = eoms( x0, 0.5*dt,    k2, t0 + 0.5*dt );
    k4 = eoms( x0,     dt,    k3, t0 +     dt );

    s = s0 + dt*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6
    v = v0 + dt*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6
  
    # State vector
    x =  [s,v]
  
    # Input to the integrator function being created
    integrator_in = C.INTEGRATOR_NUM_IN * [[]]
    integrator_in[C.INTEGRATOR_T0] = [t0]
    integrator_in[C.INTEGRATOR_TF] = [tf]
    integrator_in[C.INTEGRATOR_X0] =  x0
    integrator_in[C.INTEGRATOR_P]  =  [u,k,b]
  
    # Create a dummy state derivative vector
    xp0 = C.ssym("x",len(x)).data()
    integrator_in[C.INTEGRATOR_XP0] = xp0
  
    # Create the explicit rk4 integrator
    integrator = C.SXFunction(integrator_in,[x])
    return integrator


# Construct an explicit multistep integrator from Sundials
def create_integrator_cvodes():
    # Time 
    t = C.SX("t")

    # Differential states
    s = C.SX("s"); v = C.SX("v");
    y = [s,v]
    
    # Control
    u = C.SX("u")

    # Parameters
    k = C.SX("k")
    b = C.SX("b")
    m = 1.0   # mass

    # Differential equation - explicit form
    sdot = v
    vdot = (u - k*s - b*v)/m
    rhs = [sdot,vdot]
    
    # Input of the ode rhs
    ffcn_in = C.ODE_NUM_IN * [[]]
    ffcn_in[C.ODE_T] = [t]
    ffcn_in[C.ODE_Y] =  y
    ffcn_in[C.ODE_P] = [u, k, b]
    
    # ODE right hand side
    ffcn = C.SXFunction(ffcn_in,[rhs])
    ffcn.setOption("name","ODE right hand side")
    
    # Explicit integrator (CVODES)
    integrator = C.CVodesIntegrator(ffcn)
    
    # Set options
    integrator.setOption("exact_jacobian",True)
    #integrator.setOption("linear_multistep_method","bdf") # adams or bdf
    #integrator.setOption("nonlinear_solver_iteration","newton") # newton or functional
    integrator.setOption("fsens_err_con",True)
    integrator.setOption("abstol",1e-9)
    integrator.setOption("reltol",1e-7)
    integrator.setOption("steps_per_checkpoint",1000)
    #integrator.setOption("fsens_all_at_once",False)

    return integrator


# Main function
def idSystem(makePlots=False):
    # parameters
    k = 5.0  # spring constant
    b = 0.4  # viscous damping
    
    # noise
    Q = N.matrix( N.diag( [5e-3, 1e-2] )**2 )
    R = N.matrix( N.diag( [0.03, 0.06] )**2 )    
    
    # observation function
    def h(x):
        yOut = x
        return yOut
    
    # Time length
    T = 20.0
    
    # Shooting length
    nu = 150 # Number of control segments
    DT = N.double(T)/nu
    
    # Create integrator
    #integrator = create_integrator_cvodes()
    integrator = create_integrator_rk4()
    
    # Initialize the integrator
    integrator.init()
    
    ############# simulation
    # initial state
    s0 = 0 # initial position
    v0 = 1 # initial speed
    
    xNext = [s0,v0]
    
    Utrue = []
    Xtrue = [[s0,v0]]
    Y     = [[s0,v0]]
    time  = [0]
    
    for j in range(nu):
        u = N.sin(N.double(j)/10.0)
    
        integrator.setInput(      j*DT, C.INTEGRATOR_T0 )
        integrator.setInput(  (j+1)*DT, C.INTEGRATOR_TF )
        integrator.setInput( [u, k, b], C.INTEGRATOR_P  )
        integrator.setInput(     xNext, C.INTEGRATOR_X0 )
        integrator.evaluate()
    
        xNext = list(integrator.getOutput())
        x     = list(integrator.getOutput())
        y     = list(integrator.getOutput())
    
        # state record
        w = N.random.multivariate_normal([0,0], Q)
        x[0] += w[0]
        x[1] += w[1]
    
        # measurement
        v = N.random.multivariate_normal( [0,0], R )
        y[0] += v[0]
        y[1] += v[1]
    
        Xtrue.append(x)
        Y.append(y)
        Utrue.append(u)
        time.append(time[-1]+DT)
    
    ############## parameter estimation #################
    # Integrate over all intervals
    
    T0 = C.MX(0) # Beginning of time interval (changed from k*DT due to probable Sundials bug)
    TF = C.MX(DT) # End of time interval (changed from (k+1)*DT due to probable Sundials bug)
    
    design_variables = C.MX('design_variables', 2 + 2*(nu+1))
    k_hat = design_variables[0]
    b_hat = design_variables[1]
    
    x_hats = []
    for j in range(nu+1):
        x_hats.append(design_variables[2+2*j:2+2*j+2])
    
    logLiklihood = C.MX(0)
    
    # logLiklihood += sensor error
    for j in range(nu+1):
        yj = C.MX(2,1)
        yj[0,0] = C.MX(Y[j][0])
        yj[1,0] = C.MX(Y[j][1])
        xj = x_hats[j]
    
        # ll = 1/2 * err^T * R^-1 * err
        err = yj-h(xj)
        ll = -C.MX(0.5)*C.prod( err.T, C.prod(C.MX(R.I), err) )
        logLiklihood += ll
    
    # logLiklihood += dynamics constraint
    xdot = C.MX('xdot', 2,1)
    for j in range(nu):
        xj = x_hats[j]
        Xnext = integrator( [T0, TF, xj, C.vertcat([Utrue[j], k_hat, b_hat]), xdot, C.MX()])
    
        err = Xnext - x_hats[j+1]
        ll = -C.MX(0.5)*C.prod( err.T, C.prod( C.MX(Q.I), err) )
        logLiklihood += ll
    
    # Objective function
    F = -logLiklihood
    
    # Terminal constraints
    G = k_hat - C.MX(20.0)
    
    # Create the NLP
    ffcn = C.MXFunction([design_variables],[F]) # objective function
    gfcn = C.MXFunction([design_variables],[G]) # constraint function
    
    # Allocate an NLP solver
    #solver = C.IpoptSolver(ffcn,gfcn)
    solver = C.IpoptSolver(ffcn)
    
    # Set options
    solver.setOption("tol",1e-10)
    solver.setOption("hessian_approximation","limited-memory");
    solver.setOption("derivative_test","first-order")

    # initialize the solver
    solver.init()
    
    # Bounds on u and initial condition
    kb_min = [1, 0.1] + [-10 for j in range(2*(nu+1))] # lower bound
    solver.setInput(kb_min, C.NLP_SOLVER_LBX)
    
    kb_max = [10, 1] + [10 for j in range(2*(nu+1))] # upper bound
    solver.setInput(kb_max, C.NLP_SOLVER_UBX)
    
    guess = []
    for y in Y:
        guess.append(y[0])
        guess.append(y[1])
    kb_sol = [5, 0.4] + guess # initial guess
    solver.setInput(kb_sol, C.NLP_SOLVER_X0)
    
    # Bounds on g
    #Gmin = Gmax = [] #[10, 0]
    #solver.setInput(Gmin,C.NLP_SOLVER_LBG)
    #solver.setInput(Gmax,C.NLP_SOLVER_UBG)
    
    # Solve the problem
    solver.solve()
    
    # Get the solution
    xopt = solver.getOutput(C.NLP_SOLVER_X)
    
    # estimated parameters:
    print ""
    print "(estimated, real) k = ("+str(k)+", "+str(xopt[0])+"),\t"+str(100.0*(k-xopt[0])/k)+" % error"
    print "(estimated, real) b = ("+str(b)+", "+str(xopt[1])+"),\t"+str(100.0*(b-xopt[1])/b)+" % error"
    
    # estimated state:
    s_est = []
    v_est = []
    for j in range(0,  nu+1 ):
        s_est.append(xopt[2+2*j])
        v_est.append(xopt[2+2*j+1])

    # make plots
    if makePlots:
        plt.figure()
        plt.clf()
    
        plt.subplot(2,1,1)
        plt.xlabel('time')
        plt.ylabel('position')
        plt.plot(time, [x[0] for x in Xtrue])
        plt.plot(time, [y[0] for y in Y], '.')
        plt.plot(time, s_est)
        plt.legend(['true','meas','est'])
    
        plt.subplot(2,1,2)
        plt.xlabel('time')
        plt.ylabel('velocity')
        plt.plot(time, [x[1] for x in Xtrue])
        plt.plot(time, [y[1] for y in Y], '.')
        plt.plot(time, v_est)
        plt.legend(['true','meas','est'])
    
    #    plt.subplot(3,1,3)
    #    plt.xlabel('time')
    #    plt.ylabel('control input')
    #    plt.plot(time[:-1], Utrue, '.')
    
        plt.show()

    return {'k':xopt[0], 'b':xopt[1]}

if __name__ == '__main__':

#    idSystem(makePlots=False)
    idSystem(makePlots=True)

#    import os
#    for j in range(1000):
#        sol = idSystem(makePlots=False)
#        os.system('echo \"'+str(sol['k'])+', '+str(sol['b'])+'\" >> solutions.txt')

