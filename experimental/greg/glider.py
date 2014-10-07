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
import numpy as np
import casadi as C
import ode
import ocpSingleShooting
import ocpMultipleShooting

# main
if __name__ == '__main__':
    # setup ode

    gliderOde = ode.Ode('glider')
    gliderOde.addStates(['x','z','vx','vz'])
    gliderOde.addActions('alphaDeg')
    gliderOde.addParams('tEnd')

    def dxdt(x, u, p, t):

        # constants
        AR = 6     # aspect ration
        Cd0 = 0.03 # parasitic drag
        m = 2.0    # mass
        rho = 1.22 # air density
        A = 1.0    # area
        g = 9.8    # acceleration due to gravity

        # eom
        alpha = u['alphaDeg']*3.14159/180.0
        CL = 2.0*3.14159*alpha
        Cd = Cd0 + 1.1*CL*CL/(3.14159*AR)

        # airspeed
        norm_v = np.sqrt( x['vx']*x['vx'] + x['vz']*x['vz'] )
        
#        Lift = 0.5*rho*A*norm_v*CL*[ vz, -vx];
#        Drag = 0.5*rho*A*norm_v*Cd*[-vx, -vz];
        cAero = 0.5*rho*A*norm_v
        Fx = cAero*( CL*x['vz'] - Cd*x['vx'])
        Fz = cAero*(-CL*x['vx'] - Cd*x['vz'])
        ax = Fx/m
        az = Fz/m + g

        xDot = {}
        xDot['x'] = x['vx']
        xDot['z'] = x['vz']
        xDot['vx'] = ax
        xDot['vz'] = az

        return xDot

    gliderOde.setDxdt(dxdt)

    gliderOcp = ocpMultipleShooting.OcpMultipleShooting(gliderOde)
#    gliderOcp = ocpSingleShooting.OcpSingleShooting(gliderOde)
    gliderOcp.discretize(20)
    gliderOcp.setTimeInterval( 0.0, gliderOcp.getParams()['tEnd'] )

    # bounds
    v0 = 25
    inc0Deg = 45
    vx0 =  v0*np.cos(inc0Deg*3.14159265/180.0)
    vz0 = -v0*np.sin(inc0Deg*3.14159265/180.0)
    gliderOcp.bound(        'x',   -1,  1e6 )
    gliderOcp.bound(        'z', -1e6,    0 )
    gliderOcp.bound(       'vx',   -1,   1.1*v0 )
    gliderOcp.bound(       'vz',  -1.1*v0,   1.1*v0 )
    gliderOcp.bound( 'alphaDeg',  -20,   20 )
    gliderOcp.bound(     'tEnd',    1,   8 )    

    gliderOcp.bound(         'x',   0,   0,       0 )
    gliderOcp.bound(         'z',   0,   0,       0 )
    gliderOcp.bound(        'vx', vx0, vx0,       0 )
    gliderOcp.bound(        'vz', vz0, vz0,       0 )
    gliderOcp.bound(         'z',   0,   0, gliderOcp.N-1 )

#    for k in range(gliderOcp.N):
#        gliderOcp._addNonlconIneq( -gliderOcp.getState(k)['z'] )

    # objective function
#    F = -gliderOcp.getState(gliderOcp.N-1)['x']  # max distance
    F = -gliderOcp.getParams()['tEnd']    # max time

    ffcn = C.SXFunction([gliderOcp.designVariables],[F]) # objective function

    
    # run sim
    tEndGuess = 2.5
    uSim = 1*np.matrix(np.ones([gliderOde._Nu(), gliderOcp.N]))
    pSim = [tEndGuess]
    x0Sim = np.matrix(np.zeros([gliderOde._Nx(), 1]))
    x0Sim[2,0] = vx0
    x0Sim[3,0] = vz0

    timeSim = np.linspace(0, tEndGuess, gliderOcp.N)
    Xguess = gliderOde.runSim(timeSim, x0Sim, uSim, pSim)

#    plt.figure()
#    plt.clf()
#
#    t = np.matrix(timeSim).transpose()
#    plt.subplot(3,1,1)
#    plt.xlabel('x')
#    plt.ylabel('altitude')
#    plt.plot( Xguess[0,:].transpose(), -Xguess[1,:].transpose(), '--.')
#
#    plt.subplot(3,1,2)
#    plt.xlabel('time')
#    plt.ylabel('vx')
#    plt.plot( t, Xguess[2,:].transpose(), '.')
#
#    plt.subplot(3,1,3)
#    plt.xlabel('time')
#    plt.ylabel('vz')
#    plt.plot( t, Xguess[3,:].transpose(), '.')
#
#    plt.show()

    # initial guesses
    for k in range(0, gliderOcp.N):
        gliderOcp.guess[gliderOcp._getIdx('x',k)]  = Xguess[0,k]
        gliderOcp.guess[gliderOcp._getIdx('z',k)]  = Xguess[1,k]
        gliderOcp.guess[gliderOcp._getIdx('vx',k)] = Xguess[2,k]
        gliderOcp.guess[gliderOcp._getIdx('vz',k)] = Xguess[3,k]
        gliderOcp.guess[gliderOcp._getIdx('alphaDeg',k)] = uSim[0,k]
    
#    gliderOcp.guess[gliderOcp._getIdx('x')] = 0
#    gliderOcp.guess[gliderOcp._getIdx('z')] = 0
#    gliderOcp.guess[gliderOcp._getIdx('vx')] = 20
#    gliderOcp.guess[gliderOcp._getIdx('vz')] = -5

    gliderOcp.guess[gliderOcp._getIdx('tEnd')] = pSim[0]


    # Create the NLP
    gfcn = C.SXFunction([gliderOcp.designVariables],[C.vertcat(gliderOcp.G)]) # constraint function
    
    # Allocate an NLP solver
    solver = C.IpoptSolver(ffcn,gfcn)
    
    # Set options
#    solver.setOption("tol",1e-10)
    solver.setOption("tol",1e-13)
    solver.setOption("hessian_approximation","limited-memory");
    solver.setOption("derivative_test","first-order")

    # initialize the solver
    solver.init()
    
    # Bounds on u and initial condition
    solver.setInput(gliderOcp.lb, C.NLP_SOLVER_LBX)
    solver.setInput(gliderOcp.ub, C.NLP_SOLVER_UBX)
    solver.setInput(gliderOcp.guess, C.NLP_SOLVER_X0)
    
    # Bounds on g
    solver.setInput(gliderOcp.Gmin, C.NLP_SOLVER_LBG)
    solver.setInput(gliderOcp.Gmax, C.NLP_SOLVER_UBG)
    
    # Solve the problem
    solver.solve()
    
    # Get the solution
    xopt = solver.getOutput(C.NLP_SOLVER_X)

    print ""
    print "optimal time: "+str(xopt[gliderOcp._getIdx('tEnd')])+" seconds"

    # make plots
    time = np.linspace(0, xopt[gliderOcp._getIdx('tEnd')], gliderOcp.N)

    plt.figure()
    plt.clf()

#    t = np.matrix(timeSim).transpose()
    t = np.linspace(0, xopt[gliderOcp._getIdx('tEnd')], gliderOcp.N)
    plt.subplot(4,1,1)
    plt.xlabel('x')
    plt.ylabel('altitude')
    plt.plot( [ xopt[gliderOcp._getIdx('x',k)] for k in range(gliderOcp.N)],
              [-xopt[gliderOcp._getIdx('z',k)] for k in range(gliderOcp.N)],
              '--.')

    plt.subplot(4,1,2)
    plt.xlabel('time')
    plt.ylabel('vx')
    plt.plot( t, [-xopt[gliderOcp._getIdx('vx',k)] for k in range(gliderOcp.N)], '.')

    plt.subplot(4,1,3)
    plt.xlabel('time')
    plt.ylabel('vz')
    plt.plot( t, [-xopt[gliderOcp._getIdx('vz',k)] for k in range(gliderOcp.N)], '.')

    plt.subplot(4,1,4)
    plt.xlabel('time')
    plt.ylabel('alphaDeg')
    plt.plot( t, [-xopt[gliderOcp._getIdx('alphaDeg',k)] for k in range(gliderOcp.N)], '.')

    plt.show()
