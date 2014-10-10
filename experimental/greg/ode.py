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
import numpy as np
import casadi as C

class Ode():
    def __init__(self, name):
        self.name = name
        self.locked = False

        self.states = []
        self.actions = []
        self.params = []

    def setDxdt(self, dxdt):
        if self.locked:
            errStr = "Ode "+self.name+" has already been assigned to an Ocp and is in read-only mode"
            raise ValueError(errStr)
        self.dxdt = dxdt

    def addStates(self, states):
        if self.locked:
            errStr = "Ode "+self.name+" has already been assigned to an Ocp and is in read-only mode"
            raise ValueError(errStr)
        if not isinstance(states, list):
            states = [states]

        for x in states:
            self._assertUniqueName(x)

        self.states += states

    def addActions(self, actions):
        if self.locked:
            errStr = "Ode "+self.name+" has already been assigned to an Ocp and is in read-only mode"
            raise ValueError(errStr)
        if not isinstance(actions, list):
            actions = [actions]

        for u in actions:
            self._assertUniqueName(u)

        self.actions += actions

    def addParams(self, params):
        if self.locked:
            errStr = "Ode "+self.name+" has already been assigned to an Ocp and is in read-only mode"
            raise ValueError(errStr)
        if not isinstance(params, list):
            params = [params]

        for p in params:
            self._assertUniqueName(p)

        self.params += params

    def _Nx(self):
        return len(self.states)

    def _Nu(self):
        return len(self.actions)

    def _Nxu(self):
        return self._Nx() + self._Nu()

    def _Np(self):
        return len(self.params)

    def _stateVecToDict(self, stateVector):
        d = {}
        for (k,x) in enumerate(self.states):
            d[x] = stateVector[k]
        return d

    def _actionVecToDict(self, actionVector):
        d = {}
        for (k,u) in enumerate(self.actions):
            d[u] = actionVector[k]
        return d

    def _paramVecToDict(self, paramVector):
        d = {}
        for (k,p) in enumerate(self.params):
            d[p] = paramVector[k]
        return d

    def _stateDictToVec(self, stateDict):
        stateList = [stateDict[s] for s in self.states]

        # concatanate and return appropriate type
        if isinstance(stateList[0], C.MX): # MX
            return C.vertcat(stateList)
        elif isinstance(stateList[0], C.SX): # SX
           return C.vertcat(stateList)
        else: # numpy array
            return np.concatenate(stateList, axis=0)

    def _dxVector_dt(self, stateVector, actionVector, paramVector, t):
        x = self._stateVecToDict(stateVector)
        u = self._actionVecToDict(actionVector)
        p = self._paramVecToDict(paramVector)

        xDot = self.dxdt(x, u, p, t)

        return self._stateDictToVec(xDot)

    def _assertUniqueName(self, xup):
        if xup in self.states+self.actions+self.params:
            errStr = "Name "+xup+" is not unique"
            raise NameError(errStr)

    # inputs all are lists or vectors
    def rk4Step(self, x0Vec, u0Vec, u1Vec, pVec, t0, t1):

        dt = t1-t0 # time step

        k1 = self._dxVector_dt( x0Vec            ,             u0Vec, pVec, t0          )
        k2 = self._dxVector_dt( x0Vec + 0.5*dt*k1, 0.5*(u0Vec+u1Vec), pVec, t0 + 0.5*dt )
        k3 = self._dxVector_dt( x0Vec + 0.5*dt*k2, 0.5*(u0Vec+u1Vec), pVec, t0 + 0.5*dt )
        k4 = self._dxVector_dt( x0Vec +     dt*k3,             u1Vec, pVec, t0 +     dt )
    
        return x0Vec + dt*(k1 + 2*k2 + 2*k3 + k4)/6

#        return x0Vec + dt*k1 # euler


    # multiple rk4 steps on a single shooting interval
    def rk4Steps(self, x0Vec, u0Vec, u1Vec, pVec, t0, t1):
        N = 1
        
        x = x0Vec

        for k in range(N):
            dt = (t1 - t0)/np.double(N)
            t0_ = t0+k*dt
            t1_ = t0+(k+1)*dt

            slider0 = (t0_-t0)/(t1 - t0)
            slider1 = (t1_-t0)/(t1 - t0)

            u0_ = u0Vec*(1.0 - slider0) + u1Vec*slider0
            u1_ = u0Vec*(1.0 - slider1) + u1Vec*slider1

            x = self.rk4Step(x, u0_, u1_, pVec, t0_, t1_)

        return x


    def runSim(self, time, x0, u, p):
        N = len(time)
        X = np.matrix(np.zeros([self._Nx(), N]))
        X[:,0] = x0

        for k in range(N-1):
            u0 = u[:,k]
            u1 = u[:,k+1]
            dt = time[k+1] - time[k]
            X[:,k+1] = self.rk4Step(X[:,k], u0, u1, p, time[k], time[k+1])

        return X
