from casadi import *

# Create an NLP instance
nlp = SymbolicNLP()

# Parse an NL-file
nlp.parseNL("/home/janderss/Desktop/cuter/csfi1.nl",{"verbose":True})

# Print the NLP
#print nlp


#print nlp.x_ub
#print nlp.lambda_init


# NLP functions
ffcn = SXFunction([nlp.x],[nlp.f])
gfcn = SXFunction([nlp.x],[nlp.g])
  
# NLP solver
nlp_solver = IpoptSolver(ffcn,gfcn)
  
# Set options
# nlp_solver.setOption("max_iter",10)
#nlp_solver.setOption("verbose",True)
# nlp_solver.setOption("linear_solver","ma57")
nlp_solver.setOption("generate_hessian",True)
# nlp_solver.setOption("hessian_approximation","limited-memory")
  
# Initialize NLP solver
nlp_solver.init()
  
# Pass the bounds and initial guess
nlp_solver.setInput(nlp.x_lb,NLP_LBX)
nlp_solver.setInput(nlp.x_ub,NLP_UBX)
nlp_solver.setInput(nlp.g_lb,NLP_LBG)
nlp_solver.setInput(nlp.g_ub,NLP_UBG)
nlp_solver.setInput(nlp.x_init,NLP_X_INIT)
  
# Solve NLP
nlp_solver.solve()

