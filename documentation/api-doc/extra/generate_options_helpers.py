from casadi import *

def addExtra(metadata):

  x=SX("x")
  f = SXFunction([x],[x**2])
  f.init()
  i = IpoptSolver(f)

  for name in i.getOptionNames():
    if name in metadata["CasADi::IpoptInternal"]["options"]:
      continue
    meta = metadata["CasADi::IpoptInternal"]["options"][name] = dict()
    meta['name'] = name
    meta['type'] = i.getOptionTypeName(name)
    meta['used'] = 'CasADi::IpoptInternal'
    meta['inherit'] = False
    meta['description'] = i.getOptionDescription(name)
    try:
      meta['default'] = i.getOptionDefault(name)
    except:
      meta['default'] = ''
      pass #too bad
    #if (len(i.getOptionAllowed(name))>1):
    #  meta['description'] += "(" + "|".join(i.getOptionAllowed(name))  + ")"

  
  x=SX("x")
  f = SXFunction([x],[x**2])
  f.init()
  try:
    i = WorhpSolver(f)
  except:
    return
    
  for name in i.getOptionNames():
    meta = metadata["CasADi::WorhpInternal"]["options"][name] = dict()
    meta['name'] = name
    meta['type'] = i.getOptionTypeName(name)
    meta['used'] = 'CasADi::WorhpInternal'
    meta['inherit'] = False
    meta['description'] = i.getOptionDescription(name)
    try:
      meta['default'] = i.getOptionDefault(name)
    except:
      meta['default'] = ''
      pass #too bad
    #if (len(i.getOptionAllowed(name))>1):
    #  meta['description'] += "(" + "|".join(i.getOptionAllowed(name))  + ")"
