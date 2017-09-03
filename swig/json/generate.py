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
try:
  from lxml import etree
except:
  raise Exception("This functional relies on lxml.\nInstall it with apt-get install python-lxml.")

import sys
import re
import itertools
import json
import time

my_module = sys.argv[2]

t0 = time.time()

removed = 0
bag = dict()
# Parse the input XML
#e = etree.parse(sys.argv[1])
e = etree.iterparse(sys.argv[1],events=('end',))
for event, elem in e:
  if not elem.tag in bag:
    bag[elem.tag] = 1
  else:
    bag[elem.tag]+=1
  if elem.tag=='top':
    r = elem
  else:
    if elem.tag in ['typescopesitem']:
      removed+=1
      elem.clear()
      while elem.getprevious() is not None:
        del elem.getparent()[0]

print(bag)
print(removed)

#r = e.getroot()

def getAttribute(e,name,default=""):
  d = e.find('attributelist/attribute[@name="' + name + '"]')
  if d is None:
    return default
  return d.attrib["value"]

def getDocstring(e):
  return getAttribute(e,"feature_docstring")

def getModule(x):
  return x.xpath('ancestor-or-self::*/module/attributelist/attribute[@name="name"]')[-1].attrib['value']

def is_internal(d,msg=None):
  fd = d.find('attributelist/attribute[@name="feature_docstring"]')
  fn = d.find('attributelist/attribute[@name="name"]')
  name = [fn if fn is None else fn.attrib['value']]
  if fd is None:
    return False
    #raise Exception("is_internal find fail",d.tag,name,msg)
  v = fd.attrib['value']
  return v.strip().startswith("[INTERNAL]")
#  return "[INTERNAL]" in v
#  return not (("_EXPORT" in v) and ("CASADI_" in v))

print("elpased", time.time()-t0)
t0 = time.time()
print("get all the enums")
# get all the enums
enums = {}
for d in r.findall('*//enum'):
  if getModule(d) != my_module: continue
  sym_name = getAttribute(d,"sym_name")
  docs = getDocstring(d)
  #if is_internal(d):
  #  continue

  assert sym_name not in enums, "overwriting an enum"
  dt = enums[sym_name] = {"sym_name": sym_name, "docs": docs,"entries":{}}
  for e in d.findall('enumitem'):
    name = getAttribute(e,"name")
    ev = getAttribute(e,"enumvalueex")
    docs = getDocstring(d)
    ev = eval(ev,dict((k,v["ev"]) for k,v in dt["entries"].items()))

    assert name not in dt["entries"], "overwriting an enum entry"
    dt["entries"][name] = {"docs": docs, "ev": ev}


# get all the classes
internalClasses = []
classes0 = {}
symnameToName = {}
for c in r.findall('*//class'):
  name = c.find('attributelist/attribute[@name="name"]').attrib["value"]
  symname = c.find('attributelist/attribute[@name="sym_name"]').attrib["value"]

  if name == "casadi::"+symname:
    pass
  else:
    assert symname not in symnameToName, "overwriting a class symname"
    symnameToName[symname] = name
    assert "casadi::"+symname not in symnameToName, "overwriting a class symname"
    symnameToName["casadi::"+symname] = name

  classModule = c.find('attributelist/attribute[@name="module"]').attrib["value"]
  if name.startswith("std::vector<"): continue
  if name.startswith("std::pair<"): continue
  if classModule is None: raise ValueError("class module information missing")
  if classModule != getModule(c): raise ValueError("class module information different than getModule()")
  if classModule != my_module: continue
  if is_internal(c,msg="class"):
    internalClasses.append(name)
    continue

  assert name not in classes0, "overwriting a class"
  classes0[name] = c

internalClasses.extend(["casadi::SXNode",
                        "casadi::CallbackCPtr",
                        "casadi::DerivativeGeneratorCPtr",
                        "casadi::CustomEvaluateCPtr"])

def splitQualifiers(x):
  for q0 in ['r.','p.','q(const).']:
    if x.startswith(q0):
      (q1,ret) = splitQualifiers( x[len(q0):] )
      return ([q0]+q1, ret)
  return ([],x)

def stripQualifiers(x):
  return splitQualifiers(x[0])[1]

def internalClass(p):
  return stripQualifiers(p) in internalClasses

classes = {}
numInternalMethods = 0
numExposedMethods = 0
numExposedConstructors = 0
numInternalConstructors = 0

def getCanonicalType(x):
  (q,v) = splitQualifiers( x )
  if v in symnameToName:
    return ''.join(q)+symnameToName[v]
  else:
    return x

def getCanonicalParams(d,debug=""):
  params = []
  if d.find('attributelist/parmlist') is not None:
    for x in d.findall('attributelist/parmlist/parm/attributelist/attribute[@name="type"]'):
       name = x.findall('../attribute[@name="name"]')[0].attrib['value']
    params.append( (getCanonicalType(x.attrib['value']),"Output" if name=="OUTPUT" else "Normal") )

  return params
print("elpased", time.time()-t0)
t0 = time.time()
print("classes0")
for name,c in classes0.items():
  docs = getDocstring(c)
  if name in classes:
    data = classes[name]
  else:
    data = classes[name] = {'methods': [],"constructors":[],"docs": docs,"bases":[]}

  for d in c.findall('cdecl'):
    dname = d.find('attributelist/attribute[@name="name"]').attrib["value"]
    module = c.find('attributelist/attribute[@name="module"]').attrib["value"]
    if module != classModule:
      raise ValueError("method module",module,"!= class module",classModule)

    if (d.find('attributelist/attribute[@name="kind"]').attrib["value"]!="function"): continue
    if (d.find('attributelist/attribute[@name="kind"]').attrib["value"]!="function"): continue
    featureIgnore = d.find('attributelist/attribute[@name="feature_ignore"]')
    if featureIgnore is not None and featureIgnore.attrib['value'] == '1':
      continue

    params = getCanonicalParams(d,debug="method")

    if any([p[0].endswith("Creator") for p in params]): continue
    if any([internalClass(p[0]) for p in params]):
      numInternalMethods += 1
      continue

    rettype = getCanonicalType( d.find('attributelist/attribute[@name="type"]').attrib["value"] )

    if internalClass(rettype):
      numInternalMethods += 1
      continue
    storage = getAttribute(d,"storage")

    access = getAttribute(d,"access")
    if access=="private": continue

    if dname == "ptr": continue # WORKAROUND FOR POTENTIAL SWIG BUG

    #if is_internal(d,msg="methods"):
    #  numInternalMethods += 1
    #  continue
    numExposedMethods += 1

    docs = getDocstring(d)

    data["methods"].append((dname,params,rettype,"Static" if storage=="static" else "Normal",docs))

  for d in itertools.chain(c.findall('constructor'),c.findall('extend/constructor')):
    params = getCanonicalParams(d,debug="constructor")
    if any([p[0].endswith("Creator") for p in params]): continue
    if any([internalClass(p) for p in params]):
      numInternalConstructors += 1
      continue

    rettype = name
    access = getAttribute(d,"access")
    if access=="private": continue

    if is_internal(d,msg="constructors"):
      numInternalConstructors += 1
      continue
    if internalClass(rettype):
      numInternalConstructors += 1
      continue
    numExposedConstructors += 1
    docs = getDocstring(d)
    data["methods"].append(("CONSTRUCTOR",params,rettype,"Constructor",docs))

  for d in c.findall('attributelist/baselist/base'):
    base = d.attrib["name"]
    #if is_internal(d,msg="bases"): continue
    data["bases"].append( getCanonicalType( base ) )
print("elpased", time.time()-t0)
t0 = time.time()
print("classes")
for n,c in classes.items():
  new_bases = []
  for base in c['bases']:
    assert not base.startswith('std::'), "baseclass rename fail"
    if base.startswith('casadi::'):
      new_bases.append( getCanonicalType( base ) )
    else:
      new_bases.append( getCanonicalType( 'casadi::'+base ) )

  classes[n]['bases'] = new_bases


# get all the functions
functions = []
numInternalFunctions = 0
for d in r.findall('*//namespace/cdecl'):
  if d.find('attributelist/attribute[@name="sym_name"]') is None: continue
  if d.find('attributelist/attribute[@name="kind"]').attrib["value"]!="function": continue

  dname = d.find('attributelist/attribute[@name="sym_name"]').attrib["value"]
  name = d.find('attributelist/attribute[@name="name"]').attrib["value"]
  if dname == "dummy": continue
  code = ""
  friendwrap = "casadi_" in name

  if friendwrap:
    dname= "casadi_" + dname

  if my_module != getModule(d): continue
  if is_internal(d,msg="functions"):
    numInternalFunctions += 1
    continue

  params = getCanonicalParams(d,debug="function")

  if any([p.endswith("Creator") for p,_ in params]): continue
  if any([internalClass(p) for p,_ in params]):
    numInternalFunctions += 1
    continue

  rettype = getCanonicalType( getAttribute(d,"type") )
  if internalClass(rettype):
    numInternalFunctions += 1
    continue
  docs = getDocstring(d)

  functions.append({'funName':dname,
                    'funParams':params,
                    'funReturn':rettype,
                    'funDocs':"",#docs,
                    'funDocslink':"",
                    'funFriendwrap':friendwrap})



treedata = {"treeClasses": [],"treeFunctions": [], "treeEnums": {}}

code = sum([list(filter(lambda i: len(i.rstrip())> 0 ,x.attrib["value"].split("\n"))) for x in r.findall("*//insert/attributelist/attribute[@name='code']")],[])
treedata['treeIncludes'] = sorted(set(map(lambda x: x.rstrip(), list(filter(lambda y: re.search("^\s*#include ",y),code)))))

def getAllMethods(name,base=None):
  if base is None:
    base = name
  if name not in classes: return []
  c = classes[name]

  ret = c["methods"]
  return ret
  if name!=base:
    #ret = [(dname,params,base if mkind=="Constructor" else rettype,mkind) for (dname,params,rettype,mkind) in ret]
    # Omit baseclass constructors
    ret = list(filter(lambda x: x[3]!="Constructor", ret))

  for b in c["bases"]:
    ret = ret + getAllMethods(b,base)

  return ret
print("elpased", time.time()-t0)
t0 = time.time()
print("classes2")
myclasses = []
for k,v in classes.items():
  methods = []
  if "methods" in v:
    for (name,pars,rettype,mkind,docs) in getAllMethods(k): # v["methods"]:
      methods.append({"methodName": name, "methodReturn": rettype, "methodParams": pars, "methodKind": mkind,"methodDocs":"","methodDocslink":""})

  treedata["treeClasses"].append({"classType": k, "classMethods": methods, "classDocs": v['docs'],"classDocslink":""})
print("elpased", time.time()-t0)
t0 = time.time()
print("functions")
treedata["treeFunctions"] = functions
print("elpased", time.time()-t0)
t0 = time.time()
print("enums")
for k,v in enums.items():
  treedata["treeEnums"][k] = {
    "enumDocs": v['docs'],
    "enumDocslink": "",
    "enumEntries": dict(
       (kk , {"enumEntryDocs": vv["docs"],"enumEntryDocslink":"","enumEntryVal": vv["ev"]})
          for kk,vv in v["entries"].items())
  }
print("%5d classes %5d functions %5d enums" % (len(treedata['treeClasses']),
                                               len(treedata['treeFunctions']),
                                               len(treedata['treeEnums'])))

#print "classes:      %5d exposed %5d internal" % (len(classes),           len(internalClasses))
#print "methods:      %5d exposed %5d internal" % (numExposedMethods,      numInternalMethods)
#print "constructors: %5d exposed %5d internal" % (numExposedConstructors, numInternalConstructors)
#print "functions:    %5d exposed %5d internal" % (len(functions),         numInternalFunctions)

treedata["treeInheritance"] = dict((k, [i for i in v["bases"]]) for k,v in classes.items())
print("elpased", time.time()-t0)
t0 = time.time()

print("dump")
with open(my_module+'.json', 'w') as outfile:
    json.dump(data, outfile, indent=True)
print("elpased", time.time()-t0)
t0 = time.time()
