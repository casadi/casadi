# Get the classification of all the cuter test problems
import csv
cuter_file = csv.reader(open('cuter_problems.csv', 'rb'), delimiter=',', quotechar='|')
  
# Type of the problem objective function. Possible values:
      #N no objective function is defined,
      #C the objective function is constant,
      #L the objective function is linear,
      #Q the objective function is quadratic,
      #S the objective function is a sum of squares, and
      #O the objective function is none of the above.
allowed_objective = ['N','C','L','Q','S','O']

# Type of constraints. Possible values:
      #U the problem is unconstrained,
      #X the problem's only constraints are fixed variables,
      #B the problem's only constraints are bounds on the variables,
      #N the problem's constraints represent the adjacency matrix of a (linear) network,
      #L the problem's constraints are linear,
      #Q the problem's constraints are quadratic, and
      #O the problem's constraints are more general than any of the above alone.
allowed_constraints = ['O']

# Smoothness of the problem. Possible values:
      #R the problem is regular, that is its first and second derivatives exist and are continuous everywhere, or
      #I the problem is irregular.
allowed_smoothness = 'R'

# Filter out the problems that we are interested in:
problems = []
for (name,classification) in cuter_file:
  # If it matches our criteria
  if classification[0] in allowed_objective and classification[1] in allowed_constraints and classification[2] in allowed_smoothness:
    problems.append(name)

print problems


