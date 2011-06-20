#! KinsolSolver
#! =====================
from casadi import *
from numpy import *
from pylab import *

#! We will investigate the working of KinsolSolver with the help of the parametrically exited Duffing equation.
#!
#$ $\ddot{u}+\dot{u}-\epsilon (2 \mu \dot{u}+\alpha u^3+2 k u \cos(\Omega t))$ with $\Omega = 2 + \epsilon \sigma$. \\
#$
#$ The first order solution is $u(t)=a \cos(\frac{1}{2} \Omega t-\frac{1}{2} \gamma)$ with the modulation equations: \\
#$ $\frac{da}{d\epsilon t} = - \left[ \mu a + \frac{1}{2} k a \sin \gamma \right]$ \\
#$ $a \frac{d\gamma}{d\epsilon t} = - \left[ -\sigma a + \frac{3}{4} \alpha a^3 + k a \cos \gamma \right]$ \\
#$
#$ We seek the stationair solution to these modulation equations.

#! Parameters
eps   = SX("eps")
mu    = SX("mu")
alpha = SX("alpha")
k     = SX("k")
sigma = SX("sigma")
params = [eps,mu,alpha,k,sigma]

#! Variables
a     = SX("a")
gamma = SX("gamma")

#! Equations
res0 = mu*a+1.0/2*k*a*sin(gamma)
res1 = -sigma * a + 3.0/4*alpha*a**3+k*a*cos(gamma)

#! Numerical values
sigma_ = 0.1
alpha_ = 0.1
k_     = 0.2
params_ = [0.1,0.1,alpha_,k_,sigma_]

#! We create a KinsolSolver instance
f=SXFunction([[a,gamma],params],[[res0,res1]])
s=KinsolSolver(f)
s.setOption("strategy","linesearch")
s.setOption("abstol",1e-14)

#$ Require $a > 0$ and $\gamma < 0$
s.setOption("constraints",[2,-2])
s.init()
s.input().set(params_)

#$ Initialize [$a$,$\gamma$] with a guess and solve
s.output().set([1,-1])
s.solve()

#! Our output is:
x_ = s.output()
print "Solution = ", x_

#! Compare with the analytic solution:
x = [sqrt(4.0/3*sigma_/alpha_),-0.5*pi]
print "Reference solution = ", x

#! We show that the residual is indeed (close to) zero
f.input(0).set(s.output())
f.input(1).set(params_)
f.evaluate()
print "residual = ", f.output()

for i in range(1):
  assert(abs(x_[i]-x[i])<1e-6)

