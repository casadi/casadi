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

# Connect to remote computer
login = input('Enter remote host within string quotes: ')

# For all the problems
for fullname in glob.glob('cuter_selected' + '/*.mod'):

  # Filename
  modfile = os.path.basename(fullname)
  print 'current file is: ' + modfile
  
  # Copy source to remote computer
  os.system('scp ' + fullname + ' ' + login + ':ampl_tmp/'+ modfile)  
  
  ## Parse the ampl file
  os.system('ssh ' + login + ' "cd ampl_tmp; ampl ' + modfile + '"')

  # Name of the results file
  nlfile = os.path.splitext(modfile)[0] + '.nl'
  print 'results file should be: ' + nlfile
  
  # Get the results back
  os.system('scp ' + login + ':ampl_tmp/'+ nlfile + ' cuter_nl2/' + nlfile) 
  
  # Remove the temporaries
  os.system('ssh ' + login + ' "cd ampl_tmp; rm ' + nlfile + ' ' + modfile + '"')
  
  