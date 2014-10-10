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
import ocp
import ode
import casadi as C

class OcpSingleShooting(ocp.Ocp):

    def setTimeInterval(self, t0, tf):
        self.t0 = t0
        self.tf = tf

        #self.ssStates = C.msym("ssStates", self.ode._Nx(), self.N)
        self.ssStates = C.ssym("ssStates", self.ode._Nx(), self.N)

        self.ssStates[:,0] = self._getState0Vec()

        # dynamics constraint
        dt = (self.tf - self.t0)/(self.N - 1)
        for k in range(self.N-1):
            x0 = self._getStateVec(k)

            u0 = self._getActionVec(k)
            u1 = self._getActionVec(k+1)

            p = self._getParamVec()
    
            t0 = k*dt
            t1 = (k+1)*dt
    
            self.ssStates[:,k+1] = self.ode.rk4Step( x0, u0, u1, p, t0, t1)

    def discretize(self, N):
        if self.ode.locked:
            errStr = "Ode "+self.name+" has already been discretized and is in read-only mode"
            raise ValueError(errStr)
        self.ode.locked = True

        self.N = N

        #self.designVariables = C.msym('designVariables', self._getBigN())
        self.designVariables = C.ssym('designVariables', self._getBigN())

        self.lb = [-1e-15 for k in range(self._getBigN())]
        self.ub = [ 1e-15 for k in range(self._getBigN())]
        self.guess = [ 0 for k in range(self._getBigN())]

    def _getBigN(self): # total number of discretized states/actions/params
        return self.ode._Nx() + self.ode._Np() + self.N*self.ode._Nu()

    def _getIdx(self, xup, timeStep=None):
        if isinstance(timeStep, list):
            raise TypeError("timeStep not supposed to be a list")

        # actions
        if xup in self.ode.actions:
            if timeStep==None:
                raise ValueError("action timeStep == None")
            return self.ode._Nu()*timeStep + self.ode.actions.index(xup)

        # state0/parameters
        elif xup in self.ode.states+self.ode.params:
            if timeStep!=None:
                raise ValueError("params don't have an associated timestep")
            return self.N*self.ode._Nu() + (self.ode.states+self.ode.params).index(xup)

    def bound(self, xup, lb, ub, timeStep=None):
        # actions
        if xup in self.ode.actions:
            if timeStep == None:
                timeStep = range(0, self.N)
            elif not isinstance(timeStep, list):
                timeStep = [timeStep]
            for ts in timeStep:
                idx = self._getIdx(xup, ts)
                self.lb[idx] = lb
                self.ub[idx] = ub

        # state0
        elif xup in self.ode.states:
            if timeStep != 0:
                errStr = "Single shooting state bound timestep must be 0 (state \""+xup+"\")"
                raise ValueError(errStr)
            idx = self._getIdx(xup)
            self.lb[idx] = lb
            self.ub[idx] = ub

        # params
        elif xup in self.ode.params:
            if timeStep != None:
                errStr = "params don't have timesteps in xfull (param \""+xup+"\")"
                raise ValueError(errStr)
            idx = self._getIdx(xup)
            self.lb[idx] = lb
            self.ub[idx] = ub

        # not states/actions/params
        else:
            errStr = "bound() cannot find \""+xup+"\" in states/actions/parameters"
            raise NameError(errStr)


    def _getStateVec(self, timeStep):
        return self.ssStates[:,timeStep]

    def _getActionVec(self, timeStep):
        if timeStep >= self.N:
            errStr = "timeStep "+str(timeStep)+" is out of range (0, "+str(self.N-1)+")"
            raise ValueError(errStr)
        return self.designVariables[ timeStep*self.ode._Nu() : (timeStep+1)*self.ode._Nu() ]

    def _getState0Vec(self):
        return self.designVariables[ self.ode._Nu()*self.N : self.ode._Nu()*self.N + self.ode._Nx() ]

    def _getParamVec(self):
        idx0 = self.ode._Nu()*self.N + self.ode._Nx()
        return self.designVariables[ idx0 : idx0 + self.ode._Np() ]

    def getParams(self):
        return self.ode._paramVecToDict(self._getParamVec())
        
    def getState(self, timeStep):
        return self.ode._stateVecToDict(self._getStateVec(timeStep))
        
    def getAction(self, timeStep):
        return self.ode._actionVecToDict(self._getActionVec(timeStep))
