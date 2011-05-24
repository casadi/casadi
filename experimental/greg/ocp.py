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
        self.designVariables = C.symbolic('designVariables', self.bigBigN)
        #self.designVariables = C.MX('designVariables', self.bigBigN)


