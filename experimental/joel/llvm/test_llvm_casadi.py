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
print F.getOutput(0)
print F.getOutput(1)

