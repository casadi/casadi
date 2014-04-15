try:
  from lxml import etree
except:
  raise Exception("This functional relies on lxml.\nInstall it with apt-get install python-lxml.")

import sys
import re
import itertools
import json

my_module = sys.argv[2]

# Parse the input XML
e = etree.parse(sys.argv[1])
r = e.getroot()

def getAttribute(e,name,default=""):
  d = e.find('attributelist/attribute[@name="' + name + '"]')
  if d is None:
    return default
  return d.attrib["value"]

def getDocstring(e):
  #return ""
  return getAttribute(e,"feature_docstring")

def getModule(x):
  return x.xpath('ancestor-or-self::*/module/attributelist/attribute[@name="name"]')[-1].attrib['value']

# get all the enums
enums = {}
for d in r.findall('*//enum'):
  if getModule(d) != my_module: continue
  sym_name = getAttribute(d,"sym_name")
  docs = getDocstring(d)
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
classes = {}
for c in r.findall('*//class'):
  name = c.find('attributelist/attribute[@name="name"]').attrib["value"]
  docs = getDocstring(c)

  classModule = c.find('attributelist/attribute[@name="module"]').attrib["value"]
  if name.startswith("std::vector<"): continue
  if name.startswith("std::pair<"): continue
  if classModule is None: raise ValueError("class module information missing")
  if classModule != getModule(c): raise ValueError("class module information different than getModule()")
  if classModule != my_module: continue

  assert name not in classes, "overwriting a class"
  data = classes[name] = {'methods': [],"constructors":[],"docs": docs}

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

    if d.find('attributelist/parmlist') is None:
      params = []
    else:
      params = [ x.attrib["value"] for x in d.findall('attributelist/parmlist/parm/attributelist/attribute[@name="type"]') ]

    rettype = d.find('attributelist/attribute[@name="type"]').attrib["value"]
    storage = getAttribute(d,"storage")

    access = getAttribute(d,"access")
    if access=="private": continue

    docs = getDocstring(d)

    data["methods"].append((dname,params,rettype,"Static" if storage=="static" else "Normal",docs))

  for d in itertools.chain(c.findall('constructor'),c.findall('extend/constructor')):
    if d.find('attributelist/parmlist') is None:
      params = []
    else:
      params = [ x.attrib["value"] for x in d.findall('attributelist/parmlist/parm/attributelist/attribute[@name="type"]') ]

    rettype = name
    access = getAttribute(d,"access")
    if access=="private": continue

    docs = getDocstring(d)
    data["methods"].append((dname,params,rettype,"Constructor",docs))

  data["bases"] = []
  for d in c.findall('attributelist/baselist/base'):
    base = d.attrib["name"]
    data["bases"].append(base)

for n,c in classes.items():
  new_bases = []
  for base in c['bases']:
    assert not base.startswith('std::'), "baseclass rename fail"
    if base.startswith('casadi::'):
      new_bases.append(base)
    else:
      new_bases.append('casadi::'+base)
  classes[n]['bases'] = new_bases


# get all the functions
functions = []
for d in r.findall('*//namespace/cdecl'):
  if d.find('attributelist/attribute[@name="sym_name"]') is None: continue
  if d.find('attributelist/attribute[@name="kind"]').attrib["value"]!="function": continue

  dname = d.find('attributelist/attribute[@name="sym_name"]').attrib["value"]
  if dname == "dummy": continue
  if my_module != getModule(d): continue

  if d.find('attributelist/parmlist') is None:
    params = []
  else:
    params = [ x.attrib["value"] for x in d.findall('attributelist/parmlist/parm/attributelist/attribute[@name="type"]') ]

  rettype = getAttribute(d,"type")
  docs = getDocstring(d)

  functions.append((dname,params,rettype,docs))



treedata = {"treeClasses": [],"treeFunctions": [], "treeEnums": {}}

code = sum([filter(lambda i: len(i.rstrip())> 0 ,x.attrib["value"].split("\n")) for x in r.findall("*//insert/attributelist/attribute[@name='code']")],[])
treedata['treeIncludes'] = sorted(set(map(lambda x: x.rstrip(), filter(lambda y: re.search("^\s*#include ",y),code))))

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
    ret = filter(lambda x: x[3]!="Constructor", ret)

  for b in c["bases"]:
    ret = ret + getAllMethods(b,base)

  return ret

myclasses = []
for k,v in classes.items():
  methods = []
  if "methods" in v:
    for (name,pars,rettype,mkind,docs) in getAllMethods(k): # v["methods"]:
      methods.append({"methodName": name, "methodReturn": rettype, "methodParams": pars, "methodKind": mkind,"methodDocs":docs,"methodDocslink":""})

  treedata["treeClasses"].append({"classType": k, "classMethods": methods, "classDocs": v["docs"],"classDocslink":""})

for (name,pars,rettype,docs) in functions:
  treedata["treeFunctions"].append({"funName": name, "funReturn": rettype, "funParams": pars, "funDocs":docs,"funDocslink":""})

for k,v in enums.items():
  treedata["treeEnums"][k] = {
    "enumDocs": v['docs'],
    "enumDocslink": "",
    "enumEntries": dict(
       (kk , {"enumEntryDocs": vv["docs"],"enumEntryDocslink":"","enumEntryVal": vv["ev"]})
          for kk,vv in v["entries"].items())
  }
print "%5d classes %5d functions %5d enums" % (len(treedata['treeClasses']),
                                               len(treedata['treeFunctions']),
                                               len(treedata['treeEnums']))


treedata["treeInheritance"] = dict((k, [i for i in v["bases"]]) for k,v in classes.items())

json.dump(treedata,file(my_module+'.json','w'),indent=True)
