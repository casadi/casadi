#
#     This file is part of CasADi.
# 
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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
import unittest
from types import *
from helpers import *

class OCPtests(casadiTestCase):
  @requires("IpoptSolver")
  def testdiscrete(self):
    self.message("Linear-quadratic problem, discrete, using IPOPT")
    # inspired by www.cs.umsl.edu/~janikow/publications/1992/GAforOpt/text.pdf
    a=1.0
    b=1.0
    q=1.0
    s=1.0
    r=1.0
    x0=100
    
    N=100
    
    X=ssym("X",N+1)
    U=ssym("U",N)
    
    V = vertcat([X,U])
    
    cost = 0
    for i in range(N):
      cost = cost + s*X[i]**2+r*U[i]**2
    cost = cost + q*X[N]**2
    
    nlp = SXFunction(nlIn(x=V),nlOut(f=cost,g=vertcat([X[0]-x0,X[1:,0]-(a*X[:N,0]+b*U)])))
    
    solver = IpoptSolver(nlp)
    solver.setOption("tol",1e-5)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",100)
    solver.setOption("print_level",0)
    solver.init()
    solver.input("lbx").set([-1000 for i in range(V.size())])
    solver.input("ubx").set([1000 for i in range(V.size())])
    solver.input("lbg").set([0 for i in range(N+1)])
    solver.input("ubg").set([0 for i in range(N+1)])
    solver.solve()
    ocp_sol=solver.output("f")[0]
    # solve the ricatti equation exactly
    K = q+0.0
    for i in range(N):
      K = s+r*a**2*K/(r+b**2*K)
    exact_sol=K * x0**2
    self.assertAlmostEqual(ocp_sol,exact_sol,10,"Linear-quadratic problem solution using IPOPT")

  @requires("IpoptSolver")
  def test_singleshooting(self):
    self.message("Single shooting")
    p0 = 0.2
    y0= 1
    yc0=dy0=0
    te=0.4

    t=ssym("t")
    q=ssym("y",2,1)
    p=ssym("p",1,1)
    # y
    # y'
    f=SXFunction(daeIn(x=q,p=p,t=t),daeOut(ode=vertcat([q[1],p[0]+q[1]**2 ])))
    f.init()
    
    integrator = CVodesIntegrator(f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    integrator.setOption("verbose",False)
    integrator.setOption("steps_per_checkpoint",10000)
    integrator.setOption("t0",0)
    integrator.setOption("tf",te)

    integrator.init()

    var = MX("var",2,1)
    par = MX("par",1,1)
    parMX= par
    
    q0   = vertcat([var[0],par])
    par  = var[1]
    qend,_,_,_ = integrator.call([q0,par])
    
    parc = MX(0)
    
    f = MXFunction([var,parMX],[qend[0]])
    f.init()
    nlp = MXFunction(nlIn(x=var),nlOut(f=-f.call([var,parc])[0]))
    nlp.init()
    solver = IpoptSolver(nlp)
    solver.setOption("tol",1e-12)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",10)
    solver.setOption("derivative_test","first-order")
    solver.setOption("print_level",0)
    solver.init()
    solver.input("lbx").set([-1, -1])
    solver.input("ubx").set([1, 0.2])
    solver.solve()
    self.assertAlmostEqual(solver.output("x")[0],1,8,"X_opt")
    self.assertAlmostEqual(solver.output("x")[1],0.2,8,"X_opt")
    
    self.assertAlmostEqual(fmax(solver.output("lam_x"),0)[0],1,8,"Cost should be linear in y0")
    self.assertAlmostEqual(fmax(solver.output("lam_x"),0)[1],(sqrt(p0)*(te*yc0**2-yc0+p0*te)*tan(arctan(yc0/sqrt(p0))+sqrt(p0)*te)+yc0**2)/(2*p0*yc0**2+2*p0**2),8,"Cost should be linear in y0")
    self.assertAlmostEqual(-solver.output("f")[0],(2*y0-log(yc0**2/p0+1))/2-log(cos(arctan(yc0/sqrt(p0))+sqrt(p0)*te)),7,"Cost")
    self.assertAlmostEqual(fmax(-solver.output("lam_x"),0)[0],0,8,"Constraint is supposed to be unactive")
    self.assertAlmostEqual(fmax(-solver.output("lam_x"),0)[1],0,8,"Constraint is supposed to be unactive")

  @requires("IpoptSolver")
  def test_singleshooting2(self):
    self.message("Single shooting 2")
    p0 = 0.2
    y0= 0.2
    yc0=dy0=0.1
    te=0.4

    t=ssym("t")
    q=ssym("y",2,1)
    p=ssym("p",1,1)
    # y
    # y'
    f=SXFunction(daeIn(x=q,p=p,t=t),daeOut(ode=vertcat([q[1],p[0]+q[1]**2 ])))
    f.init()
    
    integrator = CVodesIntegrator(f)
    integrator.setOption("reltol",1e-15)
    integrator.setOption("abstol",1e-15)
    integrator.setOption("verbose",False)
    integrator.setOption("steps_per_checkpoint",10000)
    integrator.setOption("t0",0)
    integrator.setOption("tf",te)

    integrator.init()

    var = MX("var",2,1)
    par = MX("par",1,1)
    
    q0   = vertcat([var[0],par])
    parl  = var[1]
    qend,_,_,_ = integrator.call([q0,parl])
    
    parc = MX(dy0)
    
    f = MXFunction([var,par],[qend[0]])
    f.init()
    nlp = MXFunction(nlIn(x=var),nlOut(f=-f.call([var,parc])[0],g=var[0]-var[1]))
    nlp.init()
    
    solver = IpoptSolver(nlp)
    solver.setOption("tol",1e-12)
    solver.setOption("hessian_approximation", "limited-memory")
    solver.setOption("max_iter",10)
    solver.setOption("derivative_test","first-order")
    #solver.setOption("print_level",0)
    solver.init()
    solver.input("lbx").set([-1, -1])
    solver.input("ubx").set([1, 0.2])
    solver.input("lbg").set([-1])
    solver.input("ubg").set([0])
    solver.solve()

    self.assertAlmostEqual(solver.output("x")[0],0.2,6,"X_opt")
    self.assertAlmostEqual(solver.output("x")[1],0.2,6,"X_opt")
    
    self.assertAlmostEqual(fmax(solver.output("lam_x"),0)[0],0,8,"Constraint is supposed to be unactive")
    dfdp0 = (sqrt(p0)*(te*yc0**2-yc0+p0*te)*tan(arctan(yc0/sqrt(p0))+sqrt(p0)*te)+yc0**2)/(2*p0*yc0**2+2*p0**2)
    self.assertAlmostEqual(fmax(solver.output("lam_x"),0)[1],1+dfdp0,8)
    self.assertAlmostEqual(solver.output("lam_g")[0],1,8)
    self.assertAlmostEqual(-solver.output("f")[0],(2*y0-log(yc0**2/p0+1))/2-log(cos(arctan(yc0/sqrt(p0))+sqrt(p0)*te)),7,"Cost")
    self.assertAlmostEqual(fmax(-solver.output("lam_x"),0)[0],0,8,"Constraint is supposed to be unactive")
    self.assertAlmostEqual(fmax(-solver.output("lam_x"),0)[1],0,8,"Constraint is supposed to be unactive") 
    
  def test_XML(self):
    self.message("JModelica XML parsing")
    ocp = SymbolicOCP()
    ocp.parseFMI('data/cstr.xml')
    
    self.assertEqual(ocp.t0,0)
    self.assertEqual(ocp.tf,150)
    #self.assertFalse(ocp.t0_free)
    #self.assertFalse(ocp.tf_free)
    self.assertTrue(ocp.lterm.size()==0)
    self.assertTrue(ocp.mterm.size()==1)
    m = ocp.mterm
    self.assertTrue(isinstance(m,SXMatrix))
    self.assertTrue(isinstance(ocp.t,SXMatrix))
    self.assertEquals(str(m),'cost.atTime(150)')
    print dir(ocp)
    self.assertEquals(ocp.ode.size(),3)
    self.assertEquals(len(ocp.x),3) # there are three states
    c = ocp.variable("cstr.c")
    T = ocp.variable("cstr.T")
    cost = ocp.variable("cost")
    self.assertTrue(isinstance(c,Variable))
    
    self.assertEquals(c.getName(),"cstr.c")
    self.assertEquals(T.getName(),"cstr.T")
    self.assertEquals(cost.getName(),"cost")
    self.assertEquals(c.getNominal(),1000)
       
    #print c.atTime(0)
       
    c = c.var()
    T = T.var()
    cost = cost.var()
    
    u = ocp.variable("u").var()
    self.assertEquals(ocp.path.size(),3)
    #self.assertEquals(len(ocp.cfcn_lb),3)
    #self.assertEquals(len(ocp.cfcn_ub),3)
    #self.assertTrue(ocp.cfcn[0].isEqual(T)) 
    #self.assertTrue(ocp.cfcn[1].isEqual(u)) 
    #self.assertTrue(ocp.cfcn[2].isEqual(u)) 
    #self.assertTrue(ocp.cfcn_lb[0].isMinusInf()) 
    #self.assertEquals(ocp.cfcn_lb[1].getValue(),230) 
    #self.assertTrue(ocp.cfcn_lb[2].isMinusInf()) 
    #self.assertEquals(ocp.cfcn_ub[0].getValue(),350) 
    #self.assertTrue(ocp.cfcn_ub[1].isInf())
    #self.assertEquals(ocp.cfcn_ub[2].getValue(),370) 
    print ocp.initial
    print c,T,cost
    #print c.atTime(0)
    f=SXFunction([vertcat([c,T,cost])],[ocp.initial])
    f.init()
    return 
    f.evaluate()
    self.checkarray(f.output(),matrix([-956.271065,-250.051971,0]).T,"initeq")

    
    mystates = []

  @requires("IpoptSolver")
  def testMSclass_prim(self):
    self.message("CasADi multiple shooting class")
    
    ns = 20
    nx = 3
    nu = 2
    np = 0
    nh = 0
    tf = 0.2
    
    t = ssym("t")
    x0 = ssym("x0",nx)
    p = ssym("p",nu)
    xf = x0 + p[0]
    daeres = SXFunction(daeIn(t=t, x=x0, p=p),daeOut(ode=xf))
    mayer = SXFunction([x0],[7*x0[0]])
    ms = DirectMultipleShooting(daeres,mayer)
    ms.setOption("integrator",CVodesIntegrator)
    ms.setOption("number_of_grid_points",ns)
    ms.setOption("final_time",tf)
    ms.setOption("nlp_solver",IpoptSolver)
    ms.init()
    
    for i in [OCP_LBX,OCP_UBX,OCP_X_INIT]:
      self.checkarray(ms.input(i).shape,(nx,ns+1),"shape")
      
    for i in [OCP_LBU,OCP_UBU,OCP_U_INIT]:
      self.checkarray(ms.input(i).shape,(nu,ns),"shape")
    
    for i in [OCP_LBP,OCP_UBP,OCP_P_INIT]:
      self.checkarray(ms.input(i).shape,(np,1),"shape")

    for i in [OCP_LBH,OCP_UBH]:
      self.checkarray(ms.input(i).shape,(nh,ns+1),"shape")
      
    ns = 20
    nx = 3
    nu = 2
    np = 4
    nh = 2
    tf = 0.2
    
    t = ssym("t")
    x0 = ssym("x0",nx)
    p = ssym("p",nu+np)
    xf = x0 + p[0]
    daeres = SXFunction(daeIn(t=t, x=x0, p=p),daeOut(ode=xf))
    mayer = SXFunction([x0],[7*x0[0]])
    
    t = SX("t")
    cfcn = SXFunction(daeIn(t=t,x=x0, p=p),[x0[:nh,0]])
    cfcn.init()
    
    ms = DirectMultipleShooting(daeres,mayer,cfcn)
    ms.setOption("integrator",CVodesIntegrator)
    ms.setOption("number_of_grid_points",ns)
    ms.setOption("number_of_parameters",np)
    ms.setOption("final_time",tf)
    ms.setOption("nlp_solver",IpoptSolver)
    ms.init()
    
    for i in [OCP_LBX,OCP_UBX,OCP_X_INIT]:
      self.checkarray(ms.input(i).shape,(nx,ns+1),"shape")
      
    for i in [OCP_LBU,OCP_UBU,OCP_U_INIT]:
      self.checkarray(ms.input(i).shape,(nu,ns),"shape")
    
    for i in [OCP_LBP,OCP_UBP,OCP_P_INIT]:
      self.checkarray(ms.input(i).shape,(np,1),"shape")

    for i in [OCP_LBH,OCP_UBH]:
      self.checkarray(ms.input(i).shape,(nh,ns+1),"shape")

  @requires("IpoptSolver")
  def testMSclassSimple(self):
    self.message("CasADi multiple shooting class: simple example")
    """
    The problem consists of a harmonic oscilator and a power harvester.
    
    max     int_0^T u(t) * x(t)
     u,x,y
            s.t    x'(t) = y(t)
                   y'(t) = -x(t)
                   -1  <=   u(t) <= 1
                   x(0) = 1
                   y(0) = 0
                   
    The trivial solution for u(t) is u(t)=sign(x(t))
    
    
    """
    te = 2*pi
    N = 20
    t=SX("t")
    y=ssym("y",3,1)
    p=SX("p")
    f=SXFunction(daeIn(t=t, x=y, p=p),daeOut(ode=vertcat([y[1,0],-y[0,0],p*y[0,0]])))
    f.init()
    
    # Options to be passed to the integrator
    integrator_options = {}
    integrator_options["reltol"]=1e-9
    integrator_options["abstol"]=1e-9
    integrator_options["steps_per_checkpoint"]=10000
    integrator_options["t0"]=0
    integrator_options["tf"]=te/N
    
    mayer = SXFunction([y],[-y[2]])
    mayer.init()
    
    ms = DirectMultipleShooting(f,mayer)
    ms.setOption("integrator",CVodesIntegrator)
    ms.setOption("integrator_options",integrator_options)
    ms.setOption("number_of_grid_points",N)
    ms.setOption("final_time",te)
    ms.setOption("nlp_solver",IpoptSolver)
    nlp_solver_options = {}
    nlp_solver_options["tol"] = 1e-10
    nlp_solver_options["hessian_approximation"] = "limited-memory"
    nlp_solver_options["max_iter"] = 100
    nlp_solver_options["linear_solver"] = "ma57"
    nlp_solver_options["derivative_test"] = "first-order"
    ms.setOption("nlp_solver_options",nlp_solver_options)
    ms.init()
    
    ms.input("lbx").setAll(-inf)
    ms.input("ubx").setAll(inf)
    ms.input("x_init").setAll(0)
    
    ms.input("lbu").setAll(-1)
    ms.input("ubu").setAll(1)
    ms.input("u_init").setAll(0)
    
    ms.input("lbx")[0,0] = 1
    ms.input("ubx")[0,0] = 1
   
    ms.input("lbx")[1,0] = 0
    ms.input("ubx")[1,0] = 0
 
    ms.input("lbx")[2,0] = 0
    ms.input("ubx")[2,0] = 0
    
    ms.solve()
    
    self.checkarray(sign(matrix(ms.output("x_opt"))[0,:-1]),ms.output("u_opt"),"solution")


  def testMSclassSimple2(self):
    return 
    self.message("CasADi multiple shooting class: simple example 2")
    """
    The problem consists of a harmonic oscilator and a power harvester.
    
    max     x(tend)
     u,x,a
            s.t    x'(t) = a*x(t) + u(t)
                   -1  <=   u(t) <= 1
                   -2  <=   a <= 2
                   x(0) = 0
                   
    The trivial solution for u(t) is u(t)=1, a = 1
    x(t) = 1/2 (e^(2 t)-1)
    
    
    """
    te = 1
    N = 20
    t=ssym("t")
    x=ssym("x")
    a=ssym("a")
    u=ssym("u")
    f=SXFunction(daeIn(t=t, x=x, p=vertcat([a,u])),daeOut(a*x+u))
    f.init()
    
    integrator_options = {}
    integrator_options["reltol"]=1e-9
    integrator_options["abstol"]=1e-9
    integrator_options["steps_per_checkpoint"]=10000
    integrator_options["t0"]=0
    integrator_options["tf"]=te/N
    
    mayer = SXFunction([x],[-x])
    mayer.init()
    
    ms = DirectMultipleShooting(f,mayer)
    ms.setOption("integrator",CVodesIntegrator)
    ms.setOption("integrator_options",integrator_options)
    ms.setOption("number_of_grid_points",N);
    ms.setOption("number_of_parameters",1);
    ms.setOption("final_time",te);

    ms.setOption("nlp_solver",IpoptSolver)
    nlp_solver_options = {}
    nlp_solver_options["tol"] = 1e-10
    nlp_solver_options["hessian_approximation"] = "limited-memory"
    nlp_solver_options["max_iter"] = 100
    nlp_solver_options["linear_solver"] = "ma57"
    nlp_solver_options["derivative_test"] = "first-order"
    ms.setOption("nlp_solver_options",nlp_solver_options)
    ms.init()
    
    ms.input("lbx").setAll(-inf)
    ms.input("ubx").setAll(inf)
    ms.input("x_init").setAll(0)
    
    ms.input("lbu").setAll(-1)
    ms.input("ubu").setAll(1)
    ms.input("u_init").setAll(0)
    
    ms.input("lbp").setAll(-2)
    ms.input("ubp").setAll(2)
    ms.input("p_init").setAll(0)
    
    ms.input("lbx")[0,0] = 0
    ms.input("ubx")[0,0] = 0
   
    ms.solve()
    
  
if __name__ == '__main__':
    unittest.main()

