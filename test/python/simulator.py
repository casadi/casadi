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
from numpy import *
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
    integrator = Function.integrator('integrator', 'cvodes', f, opts)
    q0   = MX.sym('q0')
    par  = MX.sym('p')
    qend = integrator({'x0':q0,'p':par})['xf']
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
    opts['t0'] = 0
    opts['grid'] = tc
    integrator = Simulator('integrator', 'cvodes', f, opts)

    solution = Function('solution', {'x0':q, 'p':p, 'xf':horzcat([q*exp(t**3/(3*p)) for t in tc])},
                        Function.integrator_in(), Function.integrator_out())
    
    for f in [integrator,solution]:
      f.setInput(0.3,'x0')
      f.setInput(0.7,'p')
    
    self.checkfunction(integrator,solution,adj=False,jacobian=False,sens_der=False,evals=False,digits=6)
    
  def test_simulator_time_offset(self):
    self.message('CVodes integration: simulator time offset')
    num=self.num
    t = n.linspace(0.7,num['tend'],100)

    opts = {}
    opts['t0'] = 0.7
    opts['reltol'] = 1e-15
    opts['abstol'] = 1e-15
    opts['fsens_err_con'] = True
    opts['grid'] = t
    integrator = Simulator('integrator', 'cvodes', self.dae, opts)

    integrator.setInput([num['q0']],0)
    integrator.setInput([num['p']],1)
    integrator.evaluate()

    tend=num['tend']
    q0=num['q0']
    p=num['p']

    self.assertAlmostEqual(integrator.getOutput()[0,-1],q0*exp((tend**3-0.7**3)/(3*p)),9,'Evaluation output mismatch')
            
if __name__ == '__main__':
    unittest.main()

