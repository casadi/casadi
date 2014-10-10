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
import re

# A set of integers to be found
integers_sought = ['Total number of variables............................:',
                   'Total number of equality constraints.................:',
                   'Total number of inequality constraints...............:',
                   'Number of Iterations....:',
                   'Number of nonzeros in equality constraint Jacobian...:',
                   'Number of nonzeros in inequality constraint Jacobian.:',
                   'Number of nonzeros in Lagrangian Hessian.............:',
                   'Number of objective function evaluations             =',
                   'Number of objective gradient evaluations             =',
                   'Number of equality constraint evaluations            =',
                   'Number of inequality constraint evaluations          =',
                   'Number of inequality constraint Jacobian evaluations =',
                   'Number of Lagrangian Hessian evaluations             =']
                   
reals_sought = ['Total CPU secs in IPOPT (w/o function evaluations)   =',
                'Total CPU secs in NLP function evaluations           =']

# A set of strings to be found
strings_sought = ['EXIT:']

# Create the modified file
summary = open('summary.csv', "w" )

# Write header
summary.write('"problem", ')
for s in integers_sought + reals_sought + strings_sought:
  summary.write('"' + s + ' (ASL)", "(CasADi)", ')
summary.write('\n')

# For all the problems
for fullname in glob.glob('cuter_nl/*.nl'):

  # Filename
  nlfile = os.path.basename(fullname)
  problem = os.path.splitext(nlfile)[0]
  print 'current file is: ' + problem
  
  # Start new line in summary
  summary.write('"' + problem + '", ')

  try:

    # Solvers
    solvers = ['asl','casadi']
    text = []

    # For both the results from CasADi and ASL
    for sol in solvers:
      
      # Open the results file
      resfile=open('cuter_results/' + problem + '_' + sol + '.out', "r" )
      text.append(resfile.read())
      resfile.close()
      
    # For all the types search
    for (ii,ss) in enumerate([integers_sought,reals_sought,strings_sought]):
      
      # For all the entries
      for s in ss:
      
        # For all solvers
        for t in text:
      
          # Find the text
          index = t.find(s)
          if index>=0:
            if ii==0 or ii==1:
              # Save the numerical value
              m = re.search(r'([-+]?\d*\.\d+|\d+)',t[index:])
              if m: summary.write(m.group(0))
            else:
              # Save the output until the end of the line 
              summary.write('"' + t[index:].split('\n')[0] + '"')

          # Next column
          summary.write(', ')
        
    # Next file
    summary.write('\n')


      
  
    print "Successfully processed " + problem
  except IOError:
    print "Failed to process " + problem


# Close the summary file
summary.close()






## Get the classification of all the cuter test problems
#import csv
#cuter_file = csv.reader(open('cuter_problems.csv', 'rb'), delimiter=',', quotechar='|')
  
## Type of the problem objective function. Possible values:
      ##N no objective function is defined,
      ##C the objective function is constant,
      ##L the objective function is linear,
      ##Q the objective function is quadratic,
      ##S the objective function is a sum of squares, and
      ##O the objective function is none of the above.
#allowed_objective = ['N','C','L','Q','S','O']

## Type of constraints. Possible values:
      ##U the problem is unconstrained,
      ##X the problem's only constraints are fixed variables,
      ##B the problem's only constraints are bounds on the variables,
      ##N the problem's constraints represent the adjacency matrix of a (linear) network,
      ##L the problem's constraints are linear,
      ##Q the problem's constraints are quadratic, and
      ##O the problem's constraints are more general than any of the above alone.
#allowed_constraints = ['O']

## Smoothness of the problem. Possible values:
      ##R the problem is regular, that is its first and second derivatives exist and are continuous everywhere, or
      ##I the problem is irregular.
#allowed_smoothness = 'R'

## Filter out the problems that we are interested in:
#problems = []
#for (name,classification) in cuter_file:
  ## If it matches our criteria
  #if classification[0] in allowed_objective and classification[1] in allowed_constraints and classification[2] in allowed_smoothness:
    #problems.append(name)

## Now copy the corresponding problems to a new directory
#import shutil
#for p in problems:
  #try:

    ## Filename
    #filename = p.lower() + ".mod"
    
    ## Open the source file
    #source=open('cuter_source/'+filename, "r" )
    
    ## Create the modified file
    #destination= open('cuter_selected/'+filename, "w" )

    ## Comment out remaining lines
    #comment_out_rest = False

    ## Copy line by line
    #for line in source:
      ## Look for solve command
      #if not comment_out_rest and line.startswith("solve"):
        #comment_out_rest = True
      
      ## Comment out if necessary
      #if comment_out_rest:
        #destination.write("#")
      
      ## Write the (rest of the) line
      #destination.write(line)
      
    ## Specify filename to save to
    #destination.write("\nwrite g" + p.lower() + ";\n")
    
    ## Close the files
    #source.close()
    #destination.close()
  
    #print "Successfully copied " + filename
  #except IOError:
    #print "Failed to copy " + filename


