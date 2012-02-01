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
  os.system('scp ' + login + ':ampl_tmp/'+ nlfile + ' cuter_nl/' + nlfile) 
  
  # Remove the temporaries
  os.system('ssh ' + login + ' "cd ampl_tmp; rm ' + nlfile + ' ' + modfile + '"')
  
  