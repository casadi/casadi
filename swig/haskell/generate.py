

try:
  from lxml import etree
except:
  raise Exception("This functional relies on lxml.\nInstall it with apt-get install python-lxml.")

import sys
import re

# Parse the input XML
e = etree.parse(sys.argv[1])
r = e.getroot()

classes = {}
types_table = {}
symbol_table = {}
symbol_table_reverse = {}
enums = {}

tainted_types = {}

def ucfirst(s):
  if len(s)==0: return s
  return s[0].upper()+s[1:]



def update_types_table(p):
  if p not in types_table:
    types_table[p] = (tohaskelltype(p,alias=False),tohaskelltype(p,alias=True))

def tohaskelltype(s,alias=False,level=0):
  """
  std::vector<(CasADi::FX)>

  ->

  StdVec (CasadiClass FX)
  """
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
        return r

  if level==1:
    if s.startswith("std::vector<("):
      r = tohaskelltype(s[len("std::vector<("):-2],alias,level=2)
      if r is None: return None
      if alias:
        return r + "Vec"
      else:
        return "StdVec (" + r + ")"
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
      return "StdVec (" + r + ")"
  elif s == "std::string":
    return "StdString"
  elif s == "std::ostream":
    return "StdOstream"
  elif s == "bool":
    return "CBool"
  elif s == "std::size_t":
    return "CSize"
  elif s == "int":
    return "CInt"
  elif s == "long":
    return "CLong"
  elif s == "double":
    return "CDouble"
  elif s == "unsigned char":
    return "CUChar"
  elif s == "void":
    return "CVoid"
  elif s.startswith("std::pair"):
    return None
  elif s.startswith("CasADi::"):
    if s in  symbol_table:
      sym = symbol_table[s]
    elif  s[len("CasADi::"):] in symbol_table_reverse:
      return tohaskelltype(symbol_table_reverse[s[len("CasADi::"):]],alias,level=2)
    elif s[len("CasADi::"):] in enums:
      if alias:
        return "enum"+s[len("CasADi::"):]
      else:
        return "CasadiEnum %s" % ucfirst(s[len("CasADi::"):])
    else:
      return None

    if alias:
      return sym
    else:
      return "CasadiClass " + sym
  else:
    if s in  symbol_table:
      sym = symbol_table[s]
    elif  s in symbol_table_reverse:
      return tohaskelltype(symbol_table_reverse[s],alias,level=2)
    elif s in enums:
      if alias:
        return "enum" +s
      else:
        return "CasadiEnum %s" % ucfirst(s)
    elif s=='a().q(const).char':
      return None
    else:
      raise Exception("What is '" + s + "'")

    if alias:
      return sym
    else:
      return "CasadiClass " + sym

def getAttribute(e,name,default=""):
  d = e.find('attributelist/attribute[@name="' + name + '"]')
  if d is None:
    return default
  return d.attrib["value"]

# Build up symbol_table
for c in r.findall('*//class'):
  name = getAttribute(c,"name")
  symbol_table[name] = getAttribute(c,"sym_name",default=name)

def removens(s):
  return s.replace('CasADi::','')

symbol_table_reverse = dict([(v,k) for k,v in symbol_table.items()])
symbol_table_nonamespace = dict([(removens(k),k) for k,v in symbol_table.items()])

classes = {}

for c in r.findall('*//class'):
  name = c.find('attributelist/attribute[@name="name"]').attrib["value"]
  t = classes[name] = {}
  t["haskell"] = tohaskelltype(name)



def haskellstring(s):
  return s.replace('"','\\"').replace('\\\\"','\\"').replace("\n",'\\n')

for d in r.findall('*//enum'):
  name = getAttribute(d,"name")
  sym_name = getAttribute(d,"sym_name")
  docs = getAttribute(d,"feature_docstring")
  dt = enums[sym_name] = {"sym_name": sym_name, "docs": docs,"entries":{}}
  update_types_table(sym_name)
  for e in d.findall('enumitem'):
    name = getAttribute(e,"name")
    ev = getAttribute(e,"enumvalueex")
    docs = getAttribute(d,"feature_docstring")
    ev = eval(ev,dict((k,v["ev"]) for k,v in dt["entries"].items()))

    dt["entries"][name] = {"docs": docs, "ev": ev}

for c in r.findall('*//class'):
  name = c.find('attributelist/attribute[@name="name"]').attrib["value"]
  docs = getAttribute(c,"feature_docstring")
  data = classes[name] = {'methods': [],"constructors":[],"docs": haskellstring(docs)}



  for d in c.findall('cdecl'):
     dname = d.find('attributelist/attribute[@name="name"]').attrib["value"]
     if (d.find('attributelist/attribute[@name="kind"]').attrib["value"]!="function"): continue

     if d.find('attributelist/parmlist') is None:
       params = []
     else:
       params = [ x.attrib["value"] for x in d.findall('attributelist/parmlist/parm/attributelist/attribute[@name="type"]') ]
       for p in params:
         update_types_table(p)

     rettype = d.find('attributelist/attribute[@name="type"]').attrib["value"]
     update_types_table(rettype)
     storage = getAttribute(d,"storage")
     
     access = getAttribute(d,"access")
     if access=="private": continue

     docs = getAttribute(d,"feature_docstring")

     data["methods"].append((dname,params,rettype,"Static" if storage=="static" else "Normal",docs))

  for d in c.findall('constructor'):
     dname = symbol_table[name]
     if d.find('attributelist/parmlist') is None:
       params = []
     else:
       params = [ x.attrib["value"] for x in d.findall('attributelist/parmlist/parm/attributelist/attribute[@name="type"]') ]
       for p in params:
         update_types_table(p)

     rettype = name
     update_types_table(rettype)
     access = getAttribute(d,"access")
     if access=="private": continue
     
     docs = getAttribute(d,"feature_docstring")
     data["methods"].append((dname,params,rettype,"Constructor",docs))

  data["bases"] = []
  for d in c.findall('attributelist/baselist/base'):
    base = d.attrib["name"]
    if base in symbol_table_reverse:
      base = symbol_table_reverse[base]
      data["bases"].append(base)
    elif removens(base) in symbol_table_nonamespace:
      base = symbol_table_nonamespace[removens(base)]
      data["bases"].append(base)
    else:
      raise Exception("ouch:" + base + str(symbol_table))

functions = []
for d in r.findall('*//namespace/cdecl'):
  if d.find('attributelist/attribute[@name="sym_name"]') is None: continue
  if d.find('attributelist/attribute[@name="kind"]').attrib["value"]!="function": continue

  dname = d.find('attributelist/attribute[@name="sym_name"]').attrib["value"]
  if d.find('attributelist/parmlist') is None:
    params = []
  else:
    params = [ x.attrib["value"] for x in d.findall('attributelist/parmlist/parm/attributelist/attribute[@name="type"]') ]
    for p in params:
      update_types_table(p)

  rettype = getAttribute(d,"type")
  update_types_table(rettype)
  docs = getAttribute(d,"feature_docstring")

  functions.append((dname,params,rettype,docs))



ftree  = file('CasadiTree.hs','w')
ftree.write("{-# OPTIONS_GHC -Wall #-}\n\nmodule WriteBindings.Buildbot.CasadiTree ( classes, ioschemeclasses, tools, ioschemehelpers, enums ) where\n\n\nimport WriteBindings.Types\n\n\n")


fclasses  = file('CasadiClasses.hs','w')
fclasses.write("{-# OPTIONS_GHC -Wall #-}\n\nmodule WriteBindings.Buildbot.CasadiClasses ( CasadiEnum(..), CasadiClass(..), cppTypeCasadiPrimClass,cppTypeCasadiPrimEnum,inheritance ) where\n\n")


finclude  = file('swiginclude.hpp','w')
code = sum([filter(lambda i: len(i.rstrip())> 0 ,x.attrib["value"].split("\n")) for x in r.findall("*//insert/attributelist/attribute[@name='code']")],[])
finclude.write("\n".join(sorted(set(map(lambda x: x.rstrip(), filter(lambda y: re.search("^\s*#include ",y),code))))))

ioschemeclasses = []

def getAllMethods(name,base=None):
  if base is None:
    base = name
  if name not in classes: return []
  c = classes[name]

  ret = c["methods"]
  return ret
  if name!=base:
    #ret = [(dname,params,base if mtype=="Constructor" else rettype,mtype) for (dname,params,rettype,mtype) in ret]
    # Omit baseclass constructors
    ret = filter(lambda x: x[3]!="Constructor", ret)

  for b in c["bases"]:
    ret = ret + getAllMethods(b,base)

  return ret

myclasses = []
for k,v in classes.items():
  if ("Vector" in symbol_table[k] and "IOScheme" not in symbol_table[k]) or "Pair" in symbol_table[k]: continue
  methods = []
  if "methods" in v:
    counter = {}
    for (name,pars,rettype,mtype,docs) in getAllMethods(k): # v["methods"]:
      #print ":"+name+":",pars
      params = [types_table[p][1] for p in pars]
      if rettype not in types_table: continue
      #print ":"+name+":",pars,2,params
      t = types_table[rettype][1]
      if t is None or len(filter(lambda x : x is None, params))>0: continue
      #print ":"+name+":",pars,3
      if "VoidPointer" in name or "ptr" in name or "dummy" in name: continue
      #print ":"+name+":",pars,4
      for p in params:
        tainted_types[p] = True
      tainted_types[t] = True
      #print ":"+name+":",pars,5
      if name in counter:
        counter[name]+=1
      else:
        counter[name]=0
      methods.append("""        Method (Name "%s") %s [%s] %s (Doc "%s") """ % (name+"'"*counter[name],t,",".join(params),mtype,haskellstring(docs)))

  if name.startswith("IOScheme"):
    target = ioschemeclasses
  else:
    target = myclasses

  target.append(symbol_table[k].lower())
  ftree.write("""%s :: Class
%s = Class %s methods docs
  where
    methods =
      [
%s
      ]
    docs = Doc "%s"
      \n""" % (symbol_table[k].lower(),symbol_table[k].lower(),symbol_table[k], ",\n".join(methods),haskellstring(v["docs"])))

ftree.write("classes :: [Class]\nclasses =\n  [\n  %s\n  ]\n" % ",\n  ".join(myclasses))

tools = []
ioschemehelpers = []

counter = {}
for (name,pars,rettype,docs) in functions:
  params = [types_table[p][1] for p in pars]
  if rettype not in types_table: continue
  t = types_table[rettype][1]
  if t is None or len(filter(lambda x : x is None, params))>0: continue
  if "VoidPointer" in name or "ptr" in name or "dummy" in name: continue
  if name[0].isupper(): continue
  for p in params:
    tainted_types[p] = True
  tainted_types[t] = True
  target = tools
  if "Vec" in t and (t.endswith("Sparsity") or name.endswith("In") or name.endswith("Out")):
    #name = name[:-len("Sparsity")]
    target = ioschemehelpers
  elif "Vec" in t and t.endswith("MX") :
    #name = name[:-len("MX")]
    target = ioschemehelpers
  elif "Vec" in t and t.endswith("SX"):
    #name = name[:-len("SX")]
    target = ioschemehelpers
  if name.endswith("Struct"):
    target = ioschemehelpers
  if name in counter:
    counter[name]+=1
  else:
    counter[name]=0
  target.append("""  Function (Name "%s") %s [%s] (Doc "%s") """ % (name+"'"*counter[name],t,",".join(params),haskellstring(docs)))

ftree.write("""
tools :: [Function]
tools =
  [
%s
  ]
\n""" % ",\n".join(tools))

ftree.write("""
ioschemehelpers :: [Function]
ioschemehelpers =
  [
%s
  ]
\n""" % ",\n".join(ioschemehelpers))

ftree.write("""
ioschemeclasses :: [Class]
ioschemeclasses =
  [
  %s
  ]
\n""" % ",\n  ".join(ioschemeclasses))

enumslist = []

for k,v in enums.items():
  entrieslist = ["""("%s",Doc "%s",%d)\n""" % (kk,vv["docs"],vv["ev"]) for kk,vv in v["entries"].items()]
  enumslist.append( """CEnum %s (Doc "%s") EnumInt [\n    %s    ]\n""" % (ucfirst(k),haskellstring(v['docs']),"    ,".join(entrieslist)))
ftree.write("enums :: [CEnum]\nenums = [%s  ]\n" % "  ,".join(enumslist))


aliases = dict([(a,t) for k,(t,a) in types_table.items()])
for a,t in aliases.items():
  if t is None or a is None or not a in tainted_types: continue
  ftree.write("%s :: Type\n%s = %s\n" % (a,a,t))

exportclasses = dict([(k,v) for k,v in classes.items() if ("Vector" not in symbol_table[k] or "IOScheme" in symbol_table[k]) and "Pair" not in symbol_table[k]])

fclasses.write("data CasadiClass = %s \n  deriving (Show, Eq, Ord)\n\n\n" % ( "\n  | ".join([symbol_table[k] for k,v in exportclasses.items()])))
fclasses.write("data CasadiEnum = %s \n  deriving (Show, Eq, Ord)\n\n\n" % ( "\n  | ".join(map(lambda x: ucfirst(x),enums.keys()))))

fclasses.write('cppTypeCasadiPrimClass :: CasadiClass -> String\n');
for k,v in exportclasses.items():
  fclasses.write('cppTypeCasadiPrimClass %s = "%s"\n' % (symbol_table[k],k.replace("("," ").replace(")"," ")))

fclasses.write('cppTypeCasadiPrimEnum :: CasadiEnum -> String\n');
for k,v in enums.items():
  fclasses.write('cppTypeCasadiPrimEnum %s = "CasADi::%s"\n' % (ucfirst(k),k))

fclasses.write("\ninheritance :: [(CasadiClass,[CasadiClass])]\n")
fclasses.write(
     "inheritance = [ %s ]\n\n\n\n" % (
       "\n  , ".join(
           [
            "(%s,[%s])" % (
                    symbol_table[k],
                    ",".join([symbol_table[i] for i in v["bases"]])
              ) for k,v in exportclasses.items()])))
