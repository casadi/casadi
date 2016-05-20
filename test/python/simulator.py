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
from casadi import *
import casadi as c
import numpy
import numpy as n
import unittest
from types import *
from helpers import *



scipy_available = True
try:
	import scipy.special
	from scipy.linalg import expm
except:
	scipy_available = False

class Simulatortests(casadiTestCase):


  def setUp(self):
    # Reference solution is q0 e^((t^3-t0^3)/(3 p))
    t=SX.sym('t')
    q=SX.sym('q')
    p=SX.sym('p')
    f={'t':t, 'x':q, 'p':p, 'ode':q/p*t**2}
    opts = {}
    opts['reltol'] = 1e-15
    opts['abstol'] = 1e-15
    opts['fsens_err_con'] = True
    #opts['verbose'] = True
    opts['t0'] = 0
    opts['tf'] = 2.3
    integrator = casadi.integrator('integrator', 'cvodes', f, opts)
    q0   = MX.sym('q0')
    par  = MX.sym('p')
    qend = integrator.call({'x0':q0,'p':par})['xf']
    qe=Function('qe', [q0,par],[qend])
    self.dae = f
    self.integrator = integrator
    self.qe=qe
    self.qend=qend
    self.q0=q0
    self.par=par
    self.f = f
    self.num={'tend':2.3,'q0':7.1,'p':2}
    pass

  def test_sim_full(self):
    self.message('Simulator inputs')
    num = self.num
    N = 4
    tc = n.linspace(0,num['tend'], N)

    t=SX.sym('t')
    q=SX.sym('q')
    p=SX.sym('p')

    f={'t':t, 'x':q, 'p':p, 'ode':q/p*t**2}
    opts = {}
    opts['reltol'] = 1e-15
    opts['abstol'] = 1e-15
    opts['fsens_err_con'] = True
    #opts['verbose'] = True
    opts['grid'] = tc
    opts['output_t0'] = True
    integrator = casadi.integrator('integrator', 'cvodes', f, opts)

    solution = Function('solution', {'x0':q, 'p':p, 'xf':horzcat(*[q*exp(t**3/(3*p)) for t in tc])},
                        casadi.integrator_in(), casadi.integrator_out())
    f_in = {}
    f_in['x0']=0.3
    f_in['p']=0.7

    self.checkfunction(integrator,solution,inputs=f_in,adj=False,jacobian=False,sens_der=False,evals=False,digits=6)

  def test_simulator_time_offset(self):
    self.message('CVodes integration: simulator time offset')
    num=self.num
    t = n.linspace(0.7,num['tend'],100)

    opts = {}
    opts['reltol'] = 1e-15
    opts['abstol'] = 1e-15
    opts['fsens_err_con'] = True
    opts['grid'] = t
    opts['output_t0'] = True
    integrator = casadi.integrator('integrator', 'cvodes', self.dae, opts)

    integrator_in = [0]*integrator.n_in();integrator_in[0]=[num['q0']]
    integrator_in[1]=[num['p']]
    integrator_out = integrator.call(integrator_in)

    tend=num['tend']
    q0=num['q0']
    p=num['p']

    self.assertAlmostEqual(integrator_out[0][0,-1],q0*exp((tend**3-0.7**3)/(3*p)),9,'Evaluation output mismatch')

if __name__ == '__main__':
    unittest.main()

