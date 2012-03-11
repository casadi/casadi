from casadi import *

# Construct a simple function
x1 = ssym("x1")
x2 = ssym("x2")
r1 = sin(x2)
r2 = x1+5
  
# Create function
F = SXFunction([x1,x2],[r1,r2])
F.setOption("just_in_time",True)
F.init()

# Generate C code
F.generateCode("test.c")
  
# Pass inputs
F.setInput(10,0)
F.setInput(20,1)
  
# Evaluate
F.evaluate()
  
# Print the LLVM IR
print F

# Print results
print F.output(0)
print F.output(1)

