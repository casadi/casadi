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
from casadi.tools import *
from casadi import *
import numpy as np
from pylab import *
import os
import casadi as c
import casadi
import time


casadiAvailable = False
casadiTypes = set()
try:
  import casadi as c
  casadiAvailable = True
  casadiTypes = set([type(c.SX()),type(c.SX())])
except ImportError:
  pass
  
def TRx(a):
  constr = numpy.matrix
  if casadiAvailable and type(a) in casadiTypes:
    constr = c.SX
  return  constr([[1,0,0,0],[0,cos(a),-sin(a),0],[0,sin(a),cos(a),0],[0,0,0,1]])

def TRy(a):
  constr = numpy.matrix
  if casadiAvailable and type(a) in casadiTypes:
    constr = c.SX
  return  constr([[cos(a),0,sin(a),0],[0,1,0,0],[-sin(a),0,cos(a),0],[0,0,0,1]])

def TRz(a):
  constr = numpy.matrix
  if casadiAvailable and type(a) in casadiTypes:
    constr = c.SX
  return  constr([[cos(a),-sin(a),0,0],[sin(a),cos(a),0,0],[0,0,1,0],[0,0,0,1]])

def tr(x,y,z):
  return  numpy.matrix([[1,0,0,x],[0,1,0,y],[0,0,1,z],[0,0,0,1]])
  
def scale(a):
  return  numpy.matrix([[a,0,0],[0,a,0],[0,0,a]])
  
def Tscale(a):
  return  R2T(scale(a))
  
def Tquat(q0,q1,q2,q3):
  return R2T(quat(q0,q1,q2,q3))
  
def quat(q0,q1,q2,q3):
  """
  From Jeroen's presentation. q = [e*sin(theta/2); cos(theta/2)]
  """
  constr = numpy.matrix
  types =  set([type(q) for q in [q0,q1,q2,q3]])
  #if not(types.isdisjoint(casadiTypes)):
  #  constr = c.SX

  rho = constr([[q0],[q1],[q2]])
  rho_skew = skew(rho)
  I_3 = constr([[1.0,0,0],[0,1.0,0],[0,0,1.0]])

  #A = multiply(I_3,(numpy.dot(rho.T,-rho)+q3*q3))+numpy.dot(rho,rho.T)*2.0-q3*rho_skew*2.0
  
  b = q0
  c_ = q1
  d = q2
  a = q3
  
  a2 = a**2
  b2 = b**2
  c2 = c_**2
  d2 = d**2

  am2 = -a2
  bm2 = -b2
  cm2 = -c2
  dm2 = -d2
  
  bb = 2*b
  aa = 2*a
  
  bc2 = bb*c_
  bd2 = bb*d
  ac2 = aa*c_
  ab2 = aa*b
  ad2 = aa*d
  cd2 = 2*c_*d
  
  A = constr([[a2+b2+cm2+dm2,  bc2 - ad2,  bd2  + ac2],[bc2 + ad2, a2+bm2+c2+dm2, cd2 - ab2], [ bd2 -ac2, cd2 + ab2, a2+bm2+cm2+d2]]).T

  if not(types.isdisjoint(casadiTypes)):
    constr = c.SX
  
  return constr(A.T)

def quatOld(q0,q1,q2,q3):
  """
  From Shabana AA. Dynamics of multibody systems. Cambridge Univ Pr; 2005.
  defined as [ cos(theta/2) e*sin(theta/2) ]
  """
  constr = numpy.matrix
  types =  set([type(q) for q in [q0,q1,q2,q3]])
  #if not(types.isdisjoint(casadiTypes)):
  #  constr = c.SX
    
  E  = constr([[-q1, q0, -q3, q2],[-q2, q3, q0, -q1],[-q3,-q2,q1,q0]])
  Eb = constr([[-q1, q0, q3, -q2],[-q2, -q3, q0, q1],[-q3,q2,-q1,q0]])
  
  
  if not(types.isdisjoint(casadiTypes)):
    constr = c.SX
    
  return constr(numpy.dot(E,Eb.T))

def fullR(R_0_0,R_1_0,R_2_0, R_0_1, R_1_1, R_2_1, R_0_2, R_1_2, R_2_2):
  constr = numpy.matrix
  types =  set([type(q) for q in [R_0_0,R_1_0,R_2_0, R_0_1, R_1_1, R_2_1, R_0_2, R_1_2, R_2_2]])
  if not(types.isdisjoint(casadiTypes)):
    constr = c.SX
  return constr([[R_0_0,  R_0_1,  R_0_2],[R_1_0,  R_1_1,  R_1_2 ],[R_2_0,  R_2_1,  R_2_2 ]])
  
def TfullR(R_0_0,R_1_0,R_2_0, R_0_1, R_1_1, R_2_1, R_0_2, R_1_2, R_2_2):
  return R2T(fullR(R_0_0,R_1_0,R_2_0, R_0_1, R_1_1, R_2_1, R_0_2, R_1_2, R_2_2))
  

def origin() :
  return tr(0,0,0)
  
  
def trp(T):
  return numpy.matrix(T)[:3,3]
  
def kin_inv(T):
  R=numpy.matrix(T2R(T).T)
  constr = numpy.matrix
  if type(T) in casadiTypes:
    constr = c.SX
  return constr(vstack((hstack((R,-numpy.dot(R,trp(T)))),numpy.matrix([0,0,0,1]))))


def vectorize(vec):
  """
  Make sure the result is something you can index with single index
  """
  if hasattr(vec,"shape"):
    if vec.shape[0] > 1 and vec.shape[1] > 1:
      raise Exception("vectorize: got real matrix instead of vector like thing: %s" % str(vec))
    if vec.shape[1] > 1:
      vec = vec.T
    if hasattr(vec,"tolist"):
      vec = [ i[0] for i in vec.tolist()]
  return vec

def skew(vec):
  myvec = vectorize(vec)

  x = myvec[0]
  y = myvec[1]
  z = myvec[2]

  constr = numpy.matrix
  types =  set([type(q) for q in [x,y,z]])
  if not(types.isdisjoint(casadiTypes)):
    constr = c.SX

  return constr([[0,-z,y],[z,0,-x],[-y,x,0]])
  
def invskew(S):
  return c.SX([S[2,1],S[0,2],S[1,0]])
  
def cross(a,b):
  return c.mul(skew(a),b)

def T2R(T):
  """
   Rotational part of transformation matrix 
   
  """
  return T[0:3,0:3]
  
def R2T(R):
  """
   Pack a rotational matrix in a homogenous form
   
  """
  constr = numpy.matrix
  if type(R) in casadiTypes:
    constr = c.SX
  T  = constr([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,1.0]])
  T[:3,:3] = R
  return T
  
def T2w():
  """
   skew(w_100) = T2w(T_10)
   
  """
  
def T2W(T,p,dp):
  """
   w_101 = T2W(T_10,p,dp)
   
  """
  R = T2R(T)
  dR = c.reshape(c.mul(c.jacobian(R,p),dp),(3,3))
  return invskew(c.mul(R.T,dR))

def quatDynamics(q0,q1,q2,q3):
  """
   dot(q) = quatDynamics(q)*w_101
   
  """
  B = numpy.matrix([[q3,-q2,q1],[q2,q3,-q0],[-q1,q0,q3],[-q0,-q1,-q2]])*0.5
  return B

def T2WJ(T,p):
  """
   w_101 = T2WJ(T_10,p).diff(p,t)
   
  """
  R = T2R(T)
  RT = R.T
  
  temp = []
  for i,k in [(2,1),(0,2),(1,0)]:
     #temp.append(c.mul(c.jacobian(R[:,k],p).T,R[:,i]).T)
     temp.append(c.mul(RT[i,:],c.jacobian(R[:,k],p)))

  return c.vertcat(temp)
  
class QuadcopterAnalysis:
  def __init__(self):
    self.N = 20 # Number of control intervals
    self.tau_root = [0,0.155051,0.644949,1.000000] # Choose collocation points
    self.noiserel = 0.05 # The standard deviation of noise -to-nominal value ratio

class QuadcopterModel:
  def __init__(self,NR=4,debug=False,quatnorm=False):
    """
    Keyword arguments:
      NR    -- the number of rotors
      debug -- wether to print out debug info
      quatnorm -- add the quaternion norm to the DAE rhs
    """

    # ----------- system states and their derivatives ----
    pos = struct_symSX(["x","y","z"])     # rigid body centre of mass position [m]   {0}   
    v   = struct_symSX(["vx","vy","vz"])  # rigid body centre of mass position velocity [m/s] {0}

    NR = 4                               # Number of rotors
    
    states = struct_symSX([
      entry("p",struct=pos),
      entry("v",struct=v),
      entry("q",shape=4),                # quaternions  {0} -> {1}
      entry("w",shape=3),                # rigid body angular velocity w_101 [rad/s] {1}
      entry("r",shape=NR)                # spin speed of rotor, wrt to platform. [rad/s] Should be positive!
                                         # The signs are such that positive means lift generating, regardless of spin direction.
      
    ])
    
    pos, v, q, w, r = states[...]

    # ------------------------------------------------

    dist = struct_symSX([
      entry("Faer",shape=NR),             # Disturbance on aerodynamic forcing [N]
      entry("Caer",shape=NR)             # Disturbance on aerodynamic torques [Nm]
    ])


    # ----------------- Controls ---------------------
    controls = struct_symSX([
      entry("CR",shape=NR)              # [Nm]
          # Torques of the motors that drive the rotors, acting from platform on propeller
          # The torque signs are always positive when putting energy in the propellor,
          # regardless of spin direction.
          # 
    ])
    
    CR = controls["CR"]
    
    # ------------------------------------------------


    # ----------------  Temporary symbols --------------
    F = ssym("F",3)          # Forces acting on the platform in {1} [N]
    C = ssym("C",3)          # Torques acting on the platform in {1} [Nm]

    rotors_Faer = [ssym("Faer_%d" %i,3,1) for i in range(NR)] # Placeholder for aerodynamic force acting on propeller {1} [N]
    rotors_Caer = [ssym("Caer_%d" %i,3,1) for i in range(NR)] # Placeholder for aerodynamic torques acting on propeller {1} [Nm]

    # ---------------------------------------------------


    # ----------------- Parameters ---------------------
    
    rotor_model = struct_symSX([
         "c",        # c          Cord length [m]
         "R",        # R          Radius of propeller [m]
         "CL_alpha", # CL_alpha   Lift coefficient [-]
         "alpha_0",  # alpha_0
         "CD_alpha", # CD_alpha   Drag coefficient [-]
         "CD_i",     # CD_i       Induced drag coefficient [-]  
    ])
    
    p = struct_symSX([
      entry("rotors_model",repeat=NR,struct=rotor_model),    # Parameters that describe the rotor model
      entry("rotors_I",repeat=NR,shape=sp_diag(3)),  # Inertias of rotors [kg.m^2]
      entry("rotors_spin",repeat=NR),    # Direction of spin from each rotor. 1 means rotation around positive z.
      entry("rotors_p",repeat=NR,shape=3),  # position of rotors in {1} [m],
      entry("I",sym=casadi.diag(ssym("[Ix,Iy,Iz]"))), # Inertia of rigid body [kg.m^2]
      "m",       # Mass of the whole system [kg]
      "g",       # gravity [m/s^2]
      "rho",     # Air density [kg/m^3]
    ])
    
    I,m,g,rho = p[["I","m","g","rho"]]
 
    # --------------------------------------------------

   
    # ----------------- Parameters fillin's ---------------------

    p_ = p()
    p_["rotors_spin"] = [1,-1,1,-1]

    p_["rotors_model",:,{}] =  { "c": 0.01, "R" : 0.127, "CL_alpha": 6.0, "alpha_0": 0.15, "CD_alpha": 0.02, "CD_i": 0.05} # c          Cord length [m]

    p_["m"] = 0.5      # [kg]
    p_["g"] = 9.81     # [N/kg]
    p_["rho"] = 1.225  # [kg/m^3]

    L = 0.25
    
    I_max = p_["m"] * L**2 # Inertia of a point mass at a distance L
    I_ref = I_max/5   
    
    p_["I"] = casadi.diag([I_ref/2,I_ref/2,I_ref]) # [N.m^2]
    

    p_["rotors_p",0] = DMatrix([L,0,0])
    p_["rotors_p",1] = DMatrix([0,L,0])
    p_["rotors_p",2] = DMatrix([-L,0,0])
    p_["rotors_p",3] = DMatrix([0,-L,0])

    for i in range(NR):
        R_ = p_["rotors_model",i,"R"] #  Radius of propeller [m]
        m_ = 0.01 # Mass of a propeller [kg]
        I_max = m_ * R_**2 # Inertia of a point mass
        I_ref = I_max/5 
        p_["rotors_I",i] = casadi.diag([I_ref/2,I_ref/2,I_ref])

    if debug:
        print p.vecNZcat()
        
    dist_ = dist(0)
        
    # ----------------- Scaling ---------------------
    
    scaling_states   = states(1)
    scaling_controls = controls(1)
    
    scaling_states["r"] = 500
    scaling_controls["CR"] = 0.005
    
    scaling_dist = dist()
    
    scaling_dist["Faer"] = float(p_["m"]*p_["g"]/NR)
    scaling_dist["Caer"] = 0.0026

    # ----------- Frames ------------------
    T_10 = mul(tr(*pos),Tquat(*q))
    T_01 = kin_inv(T_10)
    R_10 = T2R(T_10)
    R_01 = R_10.T
    # -------------------------------------

    dstates = struct_symSX(states)
    
    dp,dv,dq,dw,dr = dstates[...]
    
    res = struct_SX(states) # DAE residual hand side
    # ----------- Dynamics of the body ----
    res["p"] = v - dp
    # Newton, but transform the force F from {1} to {0}
    res["v"] = mul(R_10,F) - m*dv
    # Kinematics of the quaterion.
    res["q"] = mul(quatDynamics(*q),w)-dq
    # This is a trick by Sebastien Gros to stabilize the quaternion evolution equation
    res["q"] += -q*(sumAll(q**2)-1)
    # Agular impulse H_1011
    H = mul(p["I"],w)    # Due to platform
    for i in range(NR):
      H+= mul(p["rotors_I",i], w + vertcat([0,0,p["rotors_spin",i]*r[i]])) # Due to rotor i

    dH = mul(jacobian(H,w),dw) + mul(jacobian(H,q),dq) + mul(jacobian(H,r),dr) + casadi.cross(w,H)

    res["w"] = C - dH

    for i in range(NR):
      res["r",i] = CR[i] + p["rotors_spin",i]*rotors_Caer[i][2] - p["rotors_I",i][2]*(dr[i]+dw[2]) # Dynamics of rotor i
    
    # ---------------------------------

    # Make a vector of f ?
    #if quatnorm:
    #    f = vertcat(f+[sumAll(q**2)-1])
    #else:
    #    f = vertcat(f)  

    # ------------ Force model ------------

    Fg = mul(R_01,vertcat([0,0,-g*m]))

    F_total = Fg + sum(rotors_Faer)    # Total force acting on the platform
    C_total = SX([0,0,0])                    # Total torque acting on the platform

    for i in range(NR):
       C_total[:2] += rotors_Caer[i][:2] # The x and y components propagate
       C_total[2] -= p["rotors_spin",i]*CR[i]         # the z compent moves through a serparate system
       C_total += casadi.cross(p["rotors_p",i],rotors_Faer[i]) # Torques due to thrust

    
    res = substitute(res,F,F_total)
    res = substitute(res,C,C_total)
    
    subs_before = []
    subs_after  = []
    
    v_global = mul(R_01,v)
    u_z = SX([0,0,1])
    
    # Now fill in the aerodynamic forces
    for i in range(NR):
        c,R,CL_alpha,alpha_0, CD_alpha, CD_i = p["rotors_model",i,...]
        #Bristeau P-jean, Martin P, Salaun E, Petit N. The role of propeller aerodynamics in the model of a quadrotor UAV. In: Proceedings of the European Control Conference 2009.; 2009:683-688.
        v_local = v_global + (casadi.cross(w,p["rotors_p",i])) # Velocity at rotor i
        rotors_Faer_physics =  (rho*c*R**3*r[i]**2*CL_alpha*(alpha_0/3.0-v_local[2]/(2.0*R*r[i]))) * u_z
        subs_before.append(rotors_Faer[i])
        subs_after.append(rotors_Faer_physics  + dist["Faer",i])
        rotors_Caer_physics = -p["rotors_spin",i]*rho*c*R**4*r[i]**2*(CD_alpha/4.0+CD_i*alpha_0**2*(alpha_0/4.0-2.0*v_local[2]/(3.0*r[i]*R))-CL_alpha*v_local[2]/(r[i]*R)*(alpha_0/3.0-v_local[2]/(2.0*r[i]*R))) * u_z
        subs_before.append(rotors_Caer[i])
        subs_after.append(rotors_Caer_physics  + dist["Caer",i])
    

    
    res = substitute(res,veccat(subs_before),veccat(subs_after))
    
    # Make an explicit ode
    rhs = - casadi.solve(jacobian(res,dstates),substitute(res,dstates,0))
    
    # --------------------------------------

    self.res_w = res
    self.res = substitute(res,dist,dist_)
    self.res_ = substitute(self.res,p,p_)
    
    resf = SXFunction([dstates, states, controls ],[self.res_])
    resf.init()
    self.resf = resf
    
    self.rhs_w = rhs
    
    self.rhs = substitute(rhs,dist,dist_)

    self.rhs_ = substitute(self.rhs,p,p_)

    t = SX("t")
    # We end up with a DAE that captures the system dynamics
    dae = SXFunction(daeIn(t=t,x=states,p=controls),daeOut(ode=self.rhs_))
    dae.init()
    
    self.dae = dae
    
    cdae = SXFunction(controldaeIn(t=t, x=states, u= controls,p=p),daeOut(ode=self.rhs))
    cdae.init()
    self.cdae = cdae

    self.states  = states
    self.dstates = dstates
    self.p = p
    self.p_ = p_
    self.controls = controls
    self.NR = NR
    self.w = dist
    self.w_ = dist_
    self.t = t
    
    self.states_  = states()
    self.dstates_ = states()
    self.controls_ = controls()
    
    self.scaling_states = scaling_states
    self.scaling_controls = scaling_controls
    self.scaling_dist = scaling_dist

model    = QuadcopterModel()
analysis = QuadcopterAnalysis()

controls = model.controls
states   = model.states
dstates  = model.dstates
par      = model.p
cdae     = model.cdae
res      = model.res
res_w    = model.res_w
rhs      = model.rhs
rhs_w    = model.rhs_w
dist     = model.w
dist_    = model.w_

scaling_states    = model.scaling_states
scaling_controls  = model.scaling_controls
scaling_dist      = model.scaling_dist

scaling_P = mul(scaling_states.cat,scaling_states.cat.T)
scaling_K = DMatrix.ones(controls.size,states.size)

# Number of control intervals
N = 4 #analysis.N

# Place were time is shifted
Ns = N/2

# Initial guess for time duration of one period
tf = 10


# get system dimensions
ns = states.shape[0]
nu = controls.shape[0]

waypoints = DMatrix([[-2,0,2],[2,0,1]]).T

# Formulate the estimation problem using a collocation approach

# Choose collocation points
tau_root = analysis.tau_root

# Degree of interpolating polynomial + 1
d = len(tau_root)

tau = SX("tau")


#Lf = ssym("L",ns*(ns+1)/2)
Lf = ssym("L",ns*ns)
L = SX(sp_dense(ns,ns),Lf.data())
L2P = SXFunction([Lf],[L])
L2P.init()

optvar = struct_symMX([
  (
    entry("X",repeat=[N,d],struct=states),
    entry("U",repeat=N,struct=controls),
    entry("K",repeat=N,shapestruct=(controls,states)),
    entry("P",repeat=N,shapestruct=(states,states)),
  ),
  entry("z",shape=(states.shape[0],states.shape[0]-1)),
  entry("umean"),
  entry("rmean"),
  "T",
  "Tw"
])

#print optvar["X",0,0,"p"]
#raise Exception("0")

par = struct_symMX([
  entry("Xref",shape=waypoints.shape),
  entry("model",struct=model.p)
])

T,Tw = optvar[["T","Tw"]]

print "number of variables", optvar.shape

# A vertical stack of d-by-1 Langrange functions that together form a basis on the collocation interval tau = [0..1]
Le = vertcat([numpy.prod([(tau - tau_root[r])/(tau_root[j]-tau_root[r]) for r in range(d) if not(r==j)]) for j in range(d)])
L = SXFunction([tau],[Le])
L.init()
dL = L.jacobian(0,0)
dL.init()

L.input().set(1)
L.evaluate()
Lend = DMatrix(L.output())  # Le at the end of the control interval

dLm = numSample1D(dL,DMatrix(tau_root).T)  # d-by-d

resf = SXFunction(customIO(t=model.t, x=states, dx= dstates, u=controls, p=model.p,w=model.w),[model.res_w])
resf.init()

def linear_combination(v,w):
  return sum([i*j for i,j in zip(v,w)])

### LPDE -- start

# Disturbances have a standard deviation of 5% with respect to their nominal value
Sigma = c.diag(vertcat([scaling_dist,scaling_states])**2)*analysis.noiserel

K = ssym("K",controls.size,states.size)

u = ssym("u",nu)
x = ssym("x",ns)
zj = [ ssym("z",ns) for i in range(d)]
z = vertcat(zj)

w = ssym("w",model.w.size)
G = zj[0]-x

delta = ssym("delta")


for j in range(1,d):

  dX = linear_combination(zj,dLm[:,j])
  
  [dyn] = resf.eval(customIO(
    t=0,
    x=zj[j]*scaling_states,
    dx=dX/delta*scaling_states,
    u=(-mul(K,x))*scaling_controls,
    p=model.p,
    w=w
  ))
  G.append(dyn)
  

F = mul(horzcat(zj),Lend)

Gf = SXFunction(customIO(x=x,z=z,w=w,p=model.p,K=K,delta=delta),[G,jacobian(G,x),jacobian(G,z),jacobian(G,model.w)])
Gf.init()

Ff = SXFunction([z],[F,jacobian(F,z)])
Ff.init()

linsol = CSparse(Gf.output(2).sparsity())
linsol.init()

A = msym("A",Gf.output(2).sparsity())
b = msym("b",Ff.output(1).sparsity())

linsys = MXFunction([A,b],[-linsol.solve(A,b,False)])
linsys.init()

### LPDE -- end

Ff = SXFunction([z],[F,jacobian(F,z)])
Ff.init()

par_ = par()
par_["Xref"] = waypoints

# Physical times at control intervals
ts = vertcat([c.linspace(0,Tw,Ns+1),c.linspace(Tw,T,N-Ns+1)[1:]])

# Local speed of time
dts = ts[1:]-ts[:-1]

# Physical times at collocation points
tsc = []

dLm = MX(dLm)

# Collect the equality constraints
coupling = []
collocation = []
dynamics_lpde = []
  
scaling_K = MX(scaling_K)
scaling_P = MX(scaling_P)
scaling_statesMX = MX(scaling_states)
scaling_controlsMX = MX(scaling_controls)
SigmaMX = MX(Sigma)

W = MX(DMatrix.zeros(model.w.size))

# P propagation function -- start

K = msym("K",controls.size,states.size)
P = msym("K",states.size,states.size)
u = msym("u",nu)
x = msym("x",ns)
z = msym("z",ns*d)
w = msym("w",model.w.size)
delta = msym("delta")
p = msym("p",model.p.size)

G,Gx,Gz,Gw = Gf.call(customIO(x=x,z=z,w=w,p=p,K=K,delta=delta))
    
F, Fz = Ff.call([z])

[M] = linsys.call([Gz,Fz])

Phix = mul(M,Gx)
Phiw = mul(M,horzcat([Gw,Gx]))

Pn = mul([Phix,P,Phix.T]) + mul([Phiw,SigmaMX,Phiw.T])

Pnf = MXFunction(customIO(x=x,z=z,w=w,p=p,K=K,delta=delta,P=P),[Pn])
Pnf.init()

# P propagation function -- end
  
for k in range(N):
  if k+1 < N:   #  Constraint coupling state at end of k and start of k+1
    coupling.append(optvar["X",k+1,0]-mul(optvar["X",k,horzcat],Lend))
    
    #G,Gx,Gz,Gw = Gf.call(customIO(x=optvar["X",k,0],z=optvar["X",k,vertcat],w=W,p=par["model"],K=optvar["K",k],delta=ts[k+1]-ts[k]))
    
    #F, Fz = Ff.call([optvar["X",k,vertcat]])
    
    #[M] = linsys.call([Gz,Fz])
    
    #Phix = mul(M,Gx)
    #Phiw = mul(M,horzcat([Gw,Gx]))
    
    #Pn = mul([Phix,optvar["P",k],Phix.T]) + mul([Phiw,SigmaMX,Phiw.T])
    
    [Pn] = Pnf.call(customIO(x=optvar["X",k,0],z=optvar["X",k,vertcat],w=W,p=par["model"],K=optvar["K",k],delta=ts[k+1]-ts[k],P=optvar["P",k])) 
    
    dynamics_lpde.append(optvar["P",k+1]-Pn)

  tsc.append(ts[k])
  collocation.append([])
  
  for j in range(1,d):
    tsc.append(ts[k] + (ts[k+1]-ts[k]) * tau_root[j])  # The physical time
    
    dX = linear_combination(optvar["X",k,:],dLm[:,j])
    
    [dyn] = resf.call(customIO(
      t=tsc[-1],
      x=optvar["X",k,j]*scaling_statesMX,
      dx=dX/dts[k]*scaling_statesMX,
      u=optvar["U",k]*scaling_controlsMX,
      p=par["model"],
      w=W
    ))
    collocation[-1].append(dyn)

tsf  = MXFunction([optvar],[ts])
tsf.init()
tscf = MXFunction([optvar],[vertcat(tsc)])
tscf.init()

# Memory consumption blows up when building up g

C = sumRows(states["q"]**2) - 1
Cmx = sumRows(optvar["X",0,0,"q"]**2) - 1
J = jacobian(C,states)
J = SXFunction([states],[J])
J.init()
[J] = J.call([optvar["X",0,0]])

z = optvar["z"]

g = struct_MX([
  (
    entry("coupling",expr=coupling),
    entry("collocation",expr=collocation),
    entry("dynamics_lpde",expr=dynamics_lpde),
  ),
  entry("quatnorm",expr=Cmx), # Add the quaternion norm constraint at the start
  entry("Jz",expr=mul(J,z)),   # Null space magic by Julia and Sebastien
  entry("zz",expr=mul(z.T,z)-DMatrix.eye(states.shape[0]-1)),   # Null space magic by Julia and Sebastien
  entry("p",expr=mul(z.T,optvar["X",0,0]-optvar["X",-1,-1])),    # periodicity on all states
  entry("obstacle",expr=optvar["X",:,:,lambda x: x[0]**2+x[1]**2,"p"]),
  entry("T-Tw",expr=T-Tw),
  entry("init",expr=optvar["X",0,0,"p"]*scaling_states["p"]-par["Xref",:,0]),# Initial state
  entry("waypoint",expr=optvar["X",N/2,0,"p"]*scaling_states["p"] - par["Xref",:,1])
])

# Objective function
f = T

# Favor control actions that don't deviate to much from a common mean
f += 0.01 * sumAll( (optvar["U",horzcat].T-optvar["umean"])**2 )  # control regularisation

# Favor positions close to unit quaternion
q0 = DMatrix([0,0,0,1])
f += 0.01 * sum( [ sumAll(q - q0)**2 for q in optvar["X",:,0,"q"] ] )

# Favor small angular velocities
f += 0.01 * sumAll( (optvar["X",horzcat,:,0,"w"].T)**2 )

nl = MXFunction(nlpIn(x=optvar,p=par),nlpOut(f=f,g=g))
nl.setOption("verbose",True)
nl.init()

jacG = nl.jacobian("x","g")
jacG.init()

print nl



