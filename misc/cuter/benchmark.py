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
 