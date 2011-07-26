from casadi import * 

# Implicitly defined variables
x = SX("a") 
y = SX("b")

# Parameters
p = SX("p")

# Residual function 
impfun = SXFunction([[x,y],[p]],[[x*x - y,x + 2*y + p]])

# Implicitly defined function
fun = KinsolSolver(impfun) 
fun.setLinearSolver(CSparse(CRSSparsity()))
fun.setOption("linear_solver","user_defined")
fun.setOption("abstol",1e-10)
fun.init() 

# Give the parameter a value
fun.setInput(0.1)

# Evalute
fun.setOutput([1E-7,1E-7]) # initial guess to the implicitly defined variable ([x,y])
fun.evaluate()
print fun.output()

# Change output initial guess and evaluate again
fun.setOutput([1,1],0) # initial guess to the implicitly defined variable ([x,y])
fun.evaluate()
print fun.output()

# Give forward seeds
fun.setFwdSeed(1.0)

# Give adjoint seeds
fun.setAdjSeed([0.0,1.0])

# Evaluate with sensitivities
fun.evaluate(1,1)

# Print sensitivities
print "forward sensitivities = ", fun.fwdSens()
print "adjoint sensitivities = ", fun.adjSens()



