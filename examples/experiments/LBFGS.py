import numpy as np

class LBFGS:
    def __init__(self, n, hist):
        self.n = n
        self.hist = hist
        self.y = []
        self.s = []
    
    def force_update(self, xk, xkn, gk, gkn):
        s = xkn - xk
        y = gkn - gk
        self.y = [y] + self.y[0:self.hist-1]
        self.s = [s] + self.s[0:self.hist-1]
    
    def apply(self, q, J):
        l = len(self.y)
        assert(l > 0)
        γ = -1
        α = np.nan * np.ones((l,))
        ρ = np.nan * np.ones((l,))
        for i in range(0, l):
            ρ[i] = 1. / np.dot(self.y[i][J], self.s[i][J])
            if ρ[i] <= 0: 
                continue
            α[i] = ρ[i] * np.dot(self.s[i][J], q[J])
            q[J] -= α[i] * self.y[i][J]
            if γ < 0:
                yy = np.dot(self.y[i][J], self.y[i][J])
                γ = 1. / (ρ[i] * yy)
        if γ < 0:
            return False
        q[J] *= γ
        for i in range(l - 1, -1, -1):
            if ρ[i] <= 0: 
                continue
            β = ρ[i] * np.dot(self.y[i][J], q[J])
            q[J] += self.s[i][J] * (α[i] - β)
        return True
