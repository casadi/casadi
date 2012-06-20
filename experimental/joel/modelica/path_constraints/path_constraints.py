from casadi import *
from pymodelica import compile_jmu
from pyjmi import JMUModel
import pymodelica
import os, zipfile

# Compile the modelica source
curr_dir = os.path.dirname(os.path.abspath(__file__))
jmu_name = compile_jmu("TimedVariablesTest", curr_dir+"/path_constraints.mop",'optimica','ipopt',{'generate_xml_equations':True, 'generate_fmi_me_xml':False})
sfile = zipfile.ZipFile(curr_dir+'/TimedVariablesTest.jmu','r')
mfile = sfile.extract('modelDescription.xml','.')
os.remove(curr_dir+'/TimedVariablesTest.jmu')

# Import into CasADi
ocp = casadi.SymbolicOCP()
ocp.parseFMI('modelDescription.xml',{'sort_equations':False,'eliminate_dependent':False})
ocp.sortType(True) # temporary solution: enables the new sorting
print ocp

