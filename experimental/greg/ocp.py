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

class Ocp():

    def __init__(self, ode):

        self.ode = ode

        self.G = []
        self.Gmin = []
        self.Gmax = []

    def _addNonlconIneq(self, g):
        if g.size2() != 1:
            errStr = "invalid dimensions of g "+str(g.shape)
            raise ValueError(errStr)

        self.Gmin += g.size1()*[-1e-15]
        self.Gmax += g.size1()*[0]
        self.G.append(g)

    def _addNonlconEq(self, g):
        if g.size2() != 1:
            errStr = "invalid dimensions of g "+str(g.shape)
            raise ValueError(errStr)

        self.Gmin += g.size1()*[0]
        self.Gmax += g.size1()*[0]
        self.G.append(g)

    def addNonlcon(self, lhs, rel, rhs):
        if rel == "==":
            self._addNonlconEq(lhs - rhs)
        elif rel == "<=":
            self._addNonlconIneq(lhs - rhs)
        elif rel == ">=":
            self._addNonlconIneq(rhs - lhs)
        else:
            errStr = "invalid relation \""+rel+"\""
            raise ValueError(errStr)



        


        

class MultiStageOcp():
    def __init__(self, ssOcps):
        if not isinstance(ssOcps, list):
            ssOcps = [ssOcps]

        self.ssOcps = {}
        for ssOcp in ssOcps:
            self.ssOcps[ssOcp.ode.name] = ssOcps

        self.bigBigN = sum([ssOcp._getBigN() for ssOcp in self.ssOcps.values()])
        self.designVariables = C.ssym('designVariables', self.bigBigN)
        #self.designVariables = C.MX('designVariables', self.bigBigN)


