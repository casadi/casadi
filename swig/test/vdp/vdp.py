# -*- coding: utf-8 -*-
import os
import sys
from numpy import *
import matplotlib.pyplot as plt
import zipfile

# JModelica
from jmodelica.jmi import compile_jmu
from jmodelica.jmi import JMUModel

# CasADi
from casadi import *

curr_dir = os.path.dirname(os.path.abspath(__file__));
jmu_name = compile_jmu("VDP_pack.VDP_Opt", curr_dir+"/VDP.mop",'optimica','ipopt',{'generate_xml_equations':True, 'generate_fmi_xml':False})

if True:
  vdp = JMUModel(jmu_name)
  res = vdp.optimize()

  # Extract variable profiles
  x1=res['x1']
  x2=res['x2']
  u=res['u']
  t=res['time']
  cost=res['cost']

  # Plot
  plt.figure(1)
  plt.clf()
  plt.subplot(311)
  plt.plot(t,x1)
  plt.grid()
  plt.ylabel('x1')
        
  plt.subplot(312)
  plt.plot(t,x2)
  plt.grid()
  plt.ylabel('x2')
        
  plt.subplot(313)
  plt.plot(t,u)
  plt.grid()
  plt.ylabel('u')
  plt.xlabel('time')

sfile = zipfile.ZipFile(curr_dir+'/VDP_pack_VDP_Opt.jmu','r')
mfile = sfile.extract('modelDescription.xml','.')
os.remove('VDP_pack_VDP_Opt.jmu')
os.rename('modelDescription.xml','vdp.xml')

# Allocate a parser and load the xml
parser = FMIParser('vdp.xml')

# Dump representation to screen
print "XML representation"
print parser

# Obtain the symbolic representation of the OCP
ocp = parser.parse()

# Print the ocp to screen
print ocp

# Sort the variables according to type
var = ocp.sortVariables()

# The right hand side of the ACADO functions
acado_in = ACADO_FCN_NUM_IN * [[]]

# Time
acado_in[ACADO_FCN_T] = [ocp.t]

# Convert stl vector of variables to list of expressions
def toList(v, der=False):
  ret = []
  for i in v:
    if der:
      ret.append(i.der())
    else:
      ret.append(i.sx())
  return ret

# Differential state
acado_in[ACADO_FCN_XD]  = toList(var.x)

# Algebraic state
acado_in[ACADO_FCN_XA] = toList(var.z)

# Control
acado_in[ACADO_FCN_U] = toList(var.u)

# Parameter
acado_in[ACADO_FCN_P] = toList(var.p)

# State derivative
acado_in[ACADO_FCN_XDOT] = toList(var.x,True)

# The DAE function
ffcn_out = list(ocp.dae) + list(ocp.ae)

ffcn = SXFunction(acado_in,[ffcn_out])
ffcn.setOption("ad_order",1)

# Objective function
mfcn = SXFunction(acado_in,[ocp.mterm])
mfcn.setOption("ad_order",1)

# Path constraint function
cfcn = SXFunction(acado_in,[ocp.cfcn])
cfcn.setOption("ad_order",1)
  
# Initial constraint function
rfcn = SXFunction(acado_in,[ocp.initeq])
rfcn.setOption("ad_order",1)

# Create ACADO solver
ocp_solver = AcadoInterface(ffcn,mfcn,cfcn,rfcn)

# Set options
ocp_solver.setOption("start_time",ocp.t0)
ocp_solver.setOption("final_time",ocp.tf)
num_nodes = 30
ocp_solver.setOption("number_of_shooting_nodes",num_nodes)
ocp_solver.setOption("max_num_iterations",100)
ocp_solver.setOption("kkt_tolerance",1e-4)
  
# Initialize
ocp_solver.init()

# Set bounds on states
cfcn_lb = []
for i in ocp.cfcn_lb:
  cfcn_lb.append(float(i))
ocp_solver.setInput(cfcn_lb,ACADO_LBC)
  
cfcn_ub = []
for i in ocp.cfcn_ub:
  cfcn_ub.append(float(i))
ocp_solver.setInput(cfcn_ub,ACADO_UBC)
  
# Solve the optimal control problem
ocp_solver.solve()

# Print optimal cost
cost = ocp_solver.getOutput(ACADO_COST)[0]
print "optimal cost = ", cost

# Print optimal parameters
popt = ocp_solver.getOutput(ACADO_P_OPT)
print "optimal parameter values = ", popt

# Time grid
t_opt = linspace(0,ocp.tf,num_nodes+1)

# Plot optimal control
u_opt = ocp_solver.getOutput(ACADO_U_OPT)
plt.figure(3)
plt.plot(t_opt,u_opt)

# Plot optimal state trajectory
x_opt = ocp_solver.getOutput(ACADO_X_OPT)
x_opt = array(x_opt) # create numpy array
x_opt = x_opt.reshape(num_nodes+1, 3)
plt.figure(4)
plt.plot(t_opt,x_opt)

# Show the plots
plt.ion()
plt.show()
  
