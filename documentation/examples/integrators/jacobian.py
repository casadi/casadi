#! Integrator jacobian
#! =====================
from casadi import *
from numpy import *

#! We will investigate the working of integrator jacobian with the help of the parametrically exited Duffing equation:
#!
#$ $\ddot{u}+\dot{u}-\epsilon (2 \mu \dot{u}+\alpha u^3+2 k u \cos(\Omega t))$ with $\Omega = 2 + \epsilon \sigma$.

t = SX("t")

u = SX("u") 
v = SX("v") 

eps   = SX("eps")
mu    = SX("mu")
alpha = SX("alpha")
k     = SX("k")
sigma = SX("sigma")
Omega = 2 + eps*sigma

params = [eps,mu,alpha,k,sigma]
rhs    = [v,-u-eps*(2*mu*v+alpha*u**3+2*k*u*cos(Omega*t))]

f=SXFunction({'NUM': DAE_NUM_IN, DAE_T: t, DAE_Y: [u,v], DAE_P: params},[rhs])
f.init()

integrator = CVodesIntegrator(f)

integrator.init()

#! First argument is input index, secpnd argument is output index
jac = integrator.jacobian(INTEGRATOR_P, INTEGRATOR_XF)
