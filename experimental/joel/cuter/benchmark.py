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
import os
import glob

i = 0

# Remove all results files
os.system('rm cuter_results/*.out')
os.system('rm cuter_results/*.err')

# For all the problems
for fullname in glob.glob('cuter_nl' + '/*.nl'):

  # Filename
  nlfile = os.path.basename(fullname)
  nlbase = os.path.splitext(nlfile)[0]
  print 'current file is: ' + nlfile
  
  # Solve with ASL
  os.system('timeout 300 ipopt ' + fullname + ' 1>cuter_results/'+ nlbase + '_asl.out' + ' 2>cuter_results/'+ nlbase + '_asl.err')  
  print 'ASL done'
  
  # Solve with CasADi
  os.system('timeout 300 ~/dev/casadi/trunk/build/bin/asl_reader ' + fullname + ' 1>cuter_results/'+ nlbase + '_casadi.out'+ ' 2>cuter_results/'+ nlbase + '_casadi.err')  
  print 'CasADi done'
  
  #if i>2:
    #raise Exception("a")
  #i = i+1
 