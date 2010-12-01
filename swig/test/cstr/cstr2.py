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

# Compile to XML
curr_dir = os.path.dirname(os.path.abspath(__file__));
#compiler_options = {'generate_xml_equations':True, 'generate_fmi_xml':False}
compiler_options = {'generate_xml_equations':True, 'generate_fmi_xml':False, 'enable_variable_scaling':True, 'index_reduction':True}

#jmu_name = compile_jmu("CSTRLib.Components.Two_CSTRs_stat_init", curr_dir+"/CSTRLib.mo",compiler_options=compiler_options)
jmu_name = compile_jmu("CSTR2_Opt", [curr_dir+"/CSTRLib.mo", curr_dir+"/CSTR2_Opt.mop"],compiler_options=compiler_options)

jmu_file = curr_dir + '/CSTR2_Opt.jmu'

sfile = zipfile.ZipFile(jmu_file,'r')
mfile = sfile.extract('modelDescription.xml','.')
os.remove(jmu_file)
os.rename('modelDescription.xml','cstr_init.xml')

# Parse the XML
parser = FMIParser('cstr_init.xml')
ocp = parser.parse()
ocp.sortVariables()

# Print the ocp to screen
print ocp

# Get diffential states, algebraic states and state derivatives
x = list(ocp.x)
u = list(ocp.u)
z = list(ocp.xa)[0:4]

def run_demo(with_plots=True):
    print "x = ", x
    V1 = x[0]
    CA1 = x[1]
    T1 = x[2]
    V2 = x[3]
    CA2 = x[4]
    T2 = x[5]
    u2 = x[6]
    cost = x[7]

    print "u = ", u
    der_u2c = u[0]

    print "z = ", z
    q1 = z[0]
    q2 = z[1]
    Qc = z[2]
    u1 = z[3]

    # Differential states
    der_V1  = (100-q1)
    der_CA1 = (((100/V1)-((7.2e+10*CA1)*exp(((-10000)/T1))))-((q1*CA1)/V1))-((CA1*der_V1)/V1)
    der_T1  = ((((((100*350)/V1)+(((((47800*7.2e+10)*CA1)/exp((10000/400)))*exp((((-10000)/T1)+(10000/400))))/(1000*0.239)))-((q1*T1)/V1))+Qc)-((T1*der_V1)/V1))
    der_V2  = (q1-q2)
    der_CA2 = ((((q1*CA1)/V2)-((7.2e+10*CA2)*exp(((-10000)/T2))))-((q2*CA2)/V2))-((CA2*der_V2)/V2)
    der_T2  = (((((q1*T1)/V2)+((((47800*7.2e+10)*CA2)*exp(((-10000)/T2)))/(1000*0.239)))-((q2*T2)/V2))-((T2*der_V2)/V2))
    der_u2  = der_u2c
    der_cost = (((100000*((0.03-CA1)*(0.03-CA1)))+(100000*((0.001-CA2)*(0.001-CA2))))+(50*((1-u2)*(1-u2))))

    # Algebraic states
    rhs_q1 = q1-(10*sqrt((V1-V2)))
    rhs_q2 = q2-((10*u1)*sqrt(V2))
    rhs_Qc = Qc-((-48.1909)*u2)
    rhs_u1 = u1 - 1.1

    # The ODE
    f = []
    f += [der_V1,der_CA1,der_T1,der_V2,der_CA2,der_T2,der_u2,der_cost] # Add differential equations
    f += [rhs_q1,rhs_q2,rhs_Qc,rhs_u1] # Add algebraic equations
    
    # The right hand side of the ACADO functions
    acado_in = ACADO_FCN_NUM_IN * [[]]
    
    # Time
    acado_in[ACADO_FCN_T] = list(ocp.t)

    # Differential state
    acado_in[ACADO_FCN_XD] = x

    # Algebraic state
    acado_in[ACADO_FCN_XA] = z

    # Control
    acado_in[ACADO_FCN_U] = u

    # Parameter
    acado_in[ACADO_FCN_P] = []

    print "acado_in = ", acado_in

    # The DAE function
    ffcn = SXFunction(acado_in,[f])
    ffcn.setOption("ad_order",1)

    ffcn.init()
    ffcn.setInput(0.0,ACADO_FCN_T)
    ffcn.setInput([200, 0.03, 450, 100, 0.002, 450, 1.1, 0.0],ACADO_FCN_XD)
    ffcn.setInput([100, 110, -48.1909*1.1, 1.1],ACADO_FCN_XA)
    ffcn.setInput(0.0,ACADO_FCN_U)
    ffcn.evaluate()
    print ffcn.getOutput()


    ## Objective function
    mfcn = SXFunction(acado_in,[[cost]])
    mfcn.setOption("ad_order",1)

    # Create ACADO solver
    ocp_solver = AcadoInterface(ffcn,mfcn)

    # Set options
    ocp_solver.setOption("start_time",0.0)
    ocp_solver.setOption("final_time",10.0)

    num_nodes = 10
    ocp_solver.setOption("number_of_shooting_nodes",num_nodes)
    ocp_solver.setOption("max_num_iterations",2)
    ocp_solver.setOption("max_num_integrator_steps",10000)
    ocp_solver.setOption("dynamic_sensitivity","forward_sensitivities")
    ocp_solver.setOption("kkt_tolerance",100)
    ocp_solver.setOption("absolute_tolerance",1e-4)
    ocp_solver.setOption("integrator_tolerance",1e-6)
    ocp_solver.setOption("auto_init",True)
  
    # Initialize
    ocp_solver.init()

    # Set bounds and intial guess on control
    ocp_solver.setInput([-0.1],ACADO_LBU)
    ocp_solver.setInput([ 0.1],ACADO_UBU)
    ocp_solver.setInput( (num_nodes+1)*[0],ACADO_U_GUESS)

    # Pass bounds and initial guess
    x0 = [200, 0.03, 450, 100, 0.002, 450, 1.1, 0.0]
    lb_x0 = x0 + [-inf, -inf, -inf, -inf]
    ub_x0 = x0 + [ inf,  inf,  inf,  inf]
    x_guess = x0 + [100, 100*1.1, -48.1909*1.1, 1.1]

    f += [der_V1,der_CA1,der_T1,der_V2,der_CA2,der_T2,der_u2,der_cost] # Add differential equations
    f += [rhs_q1,rhs_q2,rhs_Qc,rhs_u1] # Add algebraic equations

    ocp_solver.setInput(lb_x0,ACADO_LBX0)
    ocp_solver.setInput(ub_x0,ACADO_UBX0)
    ocp_solver.setInput((num_nodes+1) * x_guess,ACADO_X_GUESS)
  
    # Solve the optimal control problem
    ocp_solver.solve()

    # Time grid
    t_opt = linspace(0,10.0,num_nodes+1)

    # Plot optimal control
    u_opt = ocp_solver.getOutput(ACADO_U_OPT)
    u_opt = array(u_opt) # create numpy array
    u_opt = u_opt.reshape(num_nodes+1, len(u))
    plt.figure(1)
    plt.subplot(211); plt.plot(t_opt,u_opt[:,0])
    #plt.subplot(212); plt.plot(t_opt,u_opt[:,1])

    # Plot optimal state trajectory
    x_opt = ocp_solver.getOutput(ACADO_X_OPT)
    x_opt = array(x_opt) # create numpy array
    x_opt = x_opt.reshape(num_nodes+1, len(x)+len(z))
    plt.figure(2)
    plt.subplot(611); plt.plot(t_opt,x_opt[:,0])
    plt.subplot(612); plt.plot(t_opt,x_opt[:,1])
    plt.subplot(613); plt.plot(t_opt,x_opt[:,2])
    plt.subplot(614); plt.plot(t_opt,x_opt[:,3])
    plt.subplot(615); plt.plot(t_opt,x_opt[:,4])
    plt.subplot(616); plt.plot(t_opt,x_opt[:,5])

    plt.figure(3)
    plt.subplot(211); plt.plot(t_opt,x_opt[:,6])
    plt.subplot(212); plt.plot(t_opt,x_opt[:,7])



    # Show the plots
    plt.show()



    # Load a JMU model instance
    #init_model = JMUModel(jmu_name)
    
    # Set inputs for Stationary point A
    #u1_0_A = 1
    #u2_0_A = 1
    #init_model.set('u1',u1_0_A)
    #init_model.set('u2',u2_0_A)
    
    
    
    
    
    ## Solve the DAE initialization system with Ipopt
    #init_result = init_model.initialize()
    
    ## Store stationary point A
    #CA1_0_A = init_model.get('CA1')
    #CA2_0_A = init_model.get('CA2')
    #T1_0_A = init_model.get('T1')
    #T2_0_A = init_model.get('T2')
    
    ## Print some data for stationary point A
    #print(' *** Stationary point A ***')
    #print('u = [%f,%f]' % (u1_0_A,u2_0_A))
    #print('CAi = [%f,%f]' % (CA1_0_A,CA2_0_A))
    #print('Ti = [%f,%f]' % (T1_0_A,T2_0_A))
    
    ## Set inputs for stationary point B
    #u1_0_B = 1.1
    #u2_0_B = 0.9
    #init_model.set('u1',u1_0_B)
    #init_model.set('u2',u2_0_B)
    
    ## Solve the DAE initialization system with Ipopt
    #init_result = init_model.initialize()
   
    ## Stationary point B
    #CA1_0_B = init_model.get('CA1')
    #CA2_0_B = init_model.get('CA2')
    #T1_0_B = init_model.get('T1')
    #T2_0_B = init_model.get('T2')

    ## Print some data for stationary point B
    #print(' *** Stationary point B ***')
    #print('u = [%f,%f]' % (u1_0_B,u2_0_B))
    #print('CAi = [%f,%f]' % (CA1_0_B,CA2_0_B))
    #print('Ti = [%f,%f]' % (T1_0_B,T2_0_B))
    
    ### Set up and solve an optimal control problem. 

    ## Compile the Model
    #jmu_name = compile_jmu("CSTR2_Opt", 
        #[curr_dir+"/CSTRLib.mo", curr_dir+"/CSTR2_Opt.mop"],
        #compiler_options={'enable_variable_scaling':True,
            #'index_reduction':True})

    ## Load the dynamic library and XML data
    #model = JMUModel(jmu_name)

    ## Initialize the model with parameters

    ## Initialize the model to stationary point A
    #model.set('CA1_0',CA1_0_A)
    #model.set('CA2_0',CA2_0_A)
    #model.set('T1_0',T1_0_A)
    #model.set('T2_0',T2_0_A)
    
    ## Set the target values to stationary point B
    #model.set('u1_ref',u1_0_B)
    #model.set('u2_ref',u2_0_B)
    #model.set('CA1_ref',CA1_0_B)
    #model.set('CA2_ref',CA2_0_B)
    
    ## Initialize the optimization mesh
    #n_e = 50 # Number of elements 
    #hs = ones(n_e)*1./n_e # Equidistant points
    #n_cp = 3; # Number of collocation points in each element
    
    #res = model.optimize(
        #options={'n_e':n_e, 'hs':hs, 'n_cp':n_cp, 
            #'blocking_factors':2*ones(n_e/2,dtype=int), 
            #'IPOPT_options':{'tol':1e-4}})
        
    ## Extract variable profiles
    #CA1_res=res['CA1']
    #CA2_res=res['CA2']
    #T1_res=res['T1']
    #T2_res=res['T2']
    #u1_res=res['u1']
    #u2_res=res['u2']
    #der_u2_res=res['der_u2']
    
    #CA1_ref_res=res['CA1_ref']
    #CA2_ref_res=res['CA2_ref']
    
    #u1_ref_res=res['u1_ref']
    #u2_ref_res=res['u2_ref']
    
    #cost=res['cost']
    #time=res['time']

    #assert abs(cost[-1] - 1.4745648e+01) < 1e-3, \
           #"Wrong value of cost function in cstr2.py"
    
    ## Plot the results
    #if with_plots:
        #plt.figure(1)
        #plt.clf()
        #plt.hold(True)
        #plt.subplot(211)
        #plt.plot(time,CA1_res)
        #plt.plot([time[0],time[-1]],[CA1_ref_res, CA1_ref_res],'--')
        #plt.ylabel('Concentration reactor 1 [J/l]')
        #plt.grid()
        #plt.subplot(212)
        #plt.plot(time,CA2_res)
        #plt.plot([time[0],time[-1]],[CA2_ref_res, CA2_ref_res],'--')
        #plt.ylabel('Concentration reactor 2 [J/l]')
        #plt.xlabel('t [s]')
        #plt.grid()
        #plt.show()
        
        #plt.figure(2)
        #plt.clf()
        #plt.hold(True)
        #plt.subplot(211)
        #plt.plot(time,T1_res)
        #plt.ylabel('Temperature reactor 1 [K]')
        #plt.grid()
        #plt.subplot(212)
        #plt.plot(time,T2_res)
        #plt.ylabel('Temperature reactor 2 [K]')
        #plt.grid()
        #plt.xlabel('t [s]')
        #plt.show()
        
        #plt.figure(3)
        #plt.clf()
        #plt.hold(True)
        #plt.subplot(211)
        #plt.plot(time,u2_res)
        #plt.ylabel('Input 2')
        #plt.plot([time[0],time[-1]],[u2_ref_res, u2_ref_res],'--')
        #plt.grid()
        #plt.subplot(212)
        #plt.plot(time,der_u2_res)
        #plt.ylabel('Derivative of input 2')
        #plt.xlabel('t [s]')
        #plt.grid()
        #plt.show()

if __name__ == "__main__":
    run_demo()
