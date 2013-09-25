from lxml import etree
import sys

# Parse the input XML
e = etree.parse(sys.argv[1])
r = e.getroot()

classes = {}
types_table = {}
symbol_table = {}
symbol_table_reverse = {}

tainted_types = {}

def update_types_table(p):
  if p not in types_table:
    types_table[p] = (tohaskelltype(p,alias=False),tohaskelltype(p,alias=True))

def tohaskelltype(s,alias=False,level=0):
  """
  std::vector<(CasADi::FX)> 
  
  ->
  
  Vec (CasadiClass FX)
  """
  if 'IOSchemeVector' in s: return None
  if level==0:
    if s.startswith("r.q(const)."):
      r = tohaskelltype(s[len("r.q(const)."):],alias,level=1)
      if r is None: return None
      if alias:
        return "constref" + r
      else:
        return "ConstRef (" + r+ ")"
    elif s.startswith("r."):
      r = tohaskelltype(s[len("r."):],alias,level=1)
      if r is None: return None
      if alias:
        return "ref" + r
      else:
        return "Ref (" + r+ ")"
    elif s.startswith("p.q(const)."):
      r = tohaskelltype(s[len("p.q(const)."):],alias,level=1)
      if r is None: return None
      return None
      if alias:
        return "constptr" + r
      else:
        return "ConstPtr (" + r+ ")"
    elif s.startswith("q(const)."):
      return None
    elif s.startswith("p."):
      r = tohaskelltype(s[len("p."):],alias,level=1)
      if r is None: return None
      return None
      if alias:
        return "ptr" + r
      else:
        return "Ptr (" + r+ ")"
    else:
      r = tohaskelltype(s,alias,level=1)
      if r is None: return None
      if alias:
        return "val"+r
      else:
        return "Val (" + r + ")"
        
  if level==1:
    if s.startswith("std::vector<("):
      r = tohaskelltype(s[len("std::vector<("):-2],alias,level=2)
      if r is None: return None
      if alias:
        return r + "Vec"
      else:
        return "Vec (" + r + ")"
    else:
      r = tohaskelltype(s,alias,level=2)
      if r is None: return None
      if alias:
        return r
      else:
        return r
        
  if s.startswith("std::vector<("):
    r = tohaskelltype(s[len("std::vector<("):-2],alias,level=2)
    if r is None: return None
    if alias:
      return r + "Vec"
    else:
      return "Vec (" + r + ")"     
  elif s.startswith("CasADi::"):
    if s in  symbol_table:
      sym = symbol_table[s]
    elif  s[len("CasADi::"):] in symbol_table_reverse:
      return tohaskelltype(symbol_table_reverse[s[len("CasADi::"):]],alias,level=2)
    else:
      return None
    
    if alias:
      return sym
    else:
      return "NonVec (CasadiClass " + sym + ")"
  elif s == "std::string":
    return ("" if alias else "NonVec ")+ "StdString"
  elif s == "std::ostream":
    return ("" if alias else "NonVec ")+ "StdOstream"
  elif s == "bool":
    return ("" if alias else "NonVec ")+ "CBool"
  elif s == "std::size_t":
    return ("" if alias else "NonVec ")+ "CSizet"
  elif s == "int":
    return ("" if alias else "NonVec ")+ "CInt"
  elif s == "long":
    return ("" if alias else "NonVec ")+ "CLong"
  elif s == "double":
    return ("" if alias else "NonVec ")+ "CDouble"
  elif s == "unsigned char":
    return ("" if alias else "NonVec ")+ "CUChar"
  elif s == "void":
    return ("" if alias else "NonVec ")+ "CVoid"
  elif s.startswith("std::pair"):
    return None
  elif 'IOSchemeVector' in s:
    return None
  else:
    raise Exception("What is '" + s + "'")
    


# Build up symbol_table
for c in r.findall('*//class'):
  name = c.find('attributelist/attribute[@name="name"]').attrib["value"]
  if c.find('attributelist/attribute[@name="sym_name"]') is None:
    symbol_table[name] = name
  else:
    sym_name = c.find('attributelist/attribute[@name="sym_name"]').attrib["value"]
    symbol_table[name] = sym_name
    
symbol_table_reverse = dict([(v,k) for k,v in symbol_table.items()])


    
classes = {}

for c in r.findall('*//class'):
  name = c.find('attributelist/attribute[@name="name"]').attrib["value"]
  t = classes[name] = {}
  t["haskell"] = tohaskelltype(name)
  
  

for c in r.findall('*//class'):
  name = c.find('attributelist/attribute[@name="name"]').attrib["value"]
  data = classes[name] = {'methods': [],"constructors":[]}

  for d in c.findall('cdecl'):
     dname = d.find('attributelist/attribute[@name="name"]').attrib["value"]
     if d.find('attributelist/parmlist') is None:
       params = []
     else:
       params = [ x.attrib["value"] for x in d.findall('attributelist/parmlist/parm/attributelist/attribute[@name="type"]') ]
       for p in params:
         update_types_table(p)
         
     rettype = d.find('attributelist/attribute[@name="type"]').attrib["value"]
       
     data["methods"].append((dname,params,rettype,False))

  for d in c.findall('constructor'):
     dname = symbol_table[name]
     if d.find('attributelist/parmlist') is None:
       params = []
     else:
       params = [ x.attrib["value"] for x in d.findall('attributelist/parmlist/parm/attributelist/attribute[@name="type"]') ]
       for p in params:
         update_types_table(p)
         
     rettype = name
       
     data["methods"].append((dname,params,rettype,True))
    
  data["bases"] = [] 
  for d in c.findall('attributelist/baselist/base'):
    base = d.attrib["name"]
    if base in symbol_table_reverse:
      base = symbol_table_reverse[base]
      data["bases"].append(base)
    
     
functions = []
for d in r.findall('*//namespace/cdecl'):
  if d.find('attributelist/attribute[@name="sym_name"]') is None: continue
  
  dname = d.find('attributelist/attribute[@name="sym_name"]').attrib["value"]
  if d.find('attributelist/parmlist') is None:
    params = []
  else:
    params = [ x.attrib["value"] for x in d.findall('attributelist/parmlist/parm/attributelist/attribute[@name="type"]') ]
    for p in params:
      update_types_table(p)
     
  rettype = d.find('attributelist/attribute[@name="type"]').attrib["value"]
   
  functions.append((dname,params,rettype))
     

ftree  = file('CasadiTree.hs','w')
ftree.write("{-# OPTIONS_GHC -Wall #-}\n\nmodule CasadiTree ( %s, classes, tools ) where\n\n\nimport Types\n\n\n" % ",\n  ".join([symbol_table[k].lower() for k,v in classes.items() if "IOSchemeVector" not in k and "Vector" not in symbol_table[k] and "Pair" not in symbol_table[k]]))


fclasses  = file('CasadiClasses.hs','w')
fclasses.write("{-# OPTIONS_GHC -Wall#-}\n\nmodule CasadiClasses where\n\n")




def getAllMethods(name):
  if name not in classes: return []
  c = classes[name]
  
  ret = c["methods"]
  
  for b in c["bases"]:
    ret = ret + getAllMethods(b)
    
  return ret
  
myclasses = []
for k,v in classes.items():
  if "IOSchemeVector" in k: continue
  if "Vector" in symbol_table[k] or "Pair" in symbol_table[k]: continue
  methods = []
  if "methods" in v:
    counter = {}
    for (name,pars,rettype,isconstr) in getAllMethods(k): # v["methods"]:
      params = [types_table[p][1] for p in pars]
      if rettype not in types_table: continue
      t = types_table[rettype][1]
      if t is None or len(filter(lambda x : x is None, params))>0: continue
      for p in params:
        tainted_types[p] = True
      tainted_types[t] = True
      if name in counter:
        counter[name]+=1
      else:
        counter[name]=0
      methods.append("""        Method (Name "%s") %s [%s] %s""" % (name+"'"*counter[name],t,",".join(params),"Constructor" if isconstr else "Normal"))

  myclasses.append(symbol_table[k].lower())
  ftree.write("""%s :: Class
%s = Class %s methods
  where
    methods =
      [ 
%s
      ]
      \n""" % (symbol_table[k].lower(),symbol_table[k].lower(),symbol_table[k], ",\n".join(methods)))
  
ftree.write("classes :: [Class]\nclasses =\n  [\n  %s\n  ]\n" % ",\n  ".join(myclasses))
  
tools = []

counter = {}
for (name,pars,rettype) in functions:
  params = [types_table[p][1] for p in pars]
  if rettype not in types_table: continue
  t = types_table[rettype][1]
  if t is None or len(filter(lambda x : x is None, params))>0: continue
  if name[0].isupper(): continue
  for p in params:
    tainted_types[p] = True
  tainted_types[t] = True
  if name in counter:
    counter[name]+=1
  else:
    counter[name]=0
  tools.append("""  Function (Name "%s") %s [%s]""" % (name+"'"*counter[name],t,",".join(params)))
  
ftree.write("""
tools :: [Function]
tools =
  [ 
%s
  ]
\n""" % ",\n".join(tools))

  
aliases = dict([(a,t) for k,(t,a) in types_table.items()])
for a,t in aliases.items():
  if t is None or a is None or not a in tainted_types: continue
  ftree.write("%s :: Type\n%s = %s\n" % (a,a,t))
  
fclasses.write("data CasadiClass = %s \n  deriving (Show, Eq, Ord)\n\n\n" % ( "\n  | ".join([symbol_table[k] for k,v in classes.items() if "Vector" not in symbol_table[k] and "Pair" not in symbol_table[k]])))

fclasses.write('cppTypeCasadiPrim :: CasadiClass -> String\n');
for k,v in classes.items():
  if "Vector" not in symbol_table[k] and "Pair" not in symbol_table[k]:
    fclasses.write('cppTypeCasadiPrim %s = "%s"\n' % (symbol_table[k],k.replace("("," ").replace(")"," ")))

