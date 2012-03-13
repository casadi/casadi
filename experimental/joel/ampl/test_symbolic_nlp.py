from casadi import *

# Create an NLP instance
nlp = SymbolicNLP()

# Parse an NL-file
nlp.parseNL("/home/janderss/Desktop/svanberg.nl",{"verbose":True})

# Print the NLP
#print nlp


#print nlp.x_ub
#print nlp.g_ub
