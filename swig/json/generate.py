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

classes = {}
enums = {}

def ucfirst(s):
  if len(s)==0: return s
  return s[0].upper()+s[1:]

def getAttribute(e,name,default=""):
  d = e.find('attributelist/attribute[@name="' + name + '"]')
  if d is None:
    return default
  return d.attrib["value"]

classes = {}

for d in r.findall('*//enum'):
  name = getAttribute(d,"name")
  sym_name = getAttribute(d,"sym_name")
  docs = getAttribute(d,"feature_docstring")
  dt = enums[sym_name] = {"sym_name": sym_name, "docs": docs,"entries":{}}
  for e in d.findall('enumitem'):
    name = getAttribute(e,"name")
    ev = getAttribute(e,"enumvalueex")
    docs = getAttribute(d,"feature_docstring")
    ev = eval(ev,dict((k,v["ev"]) for k,v in dt["entries"].items()))

    dt["entries"][name] = {"docs": docs, "ev": ev}

for c in r.findall('*//class'):
  name = c.find('attributelist/attribute[@name="name"]').attrib["value"]
  docs = getAttribute(c,"feature_docstring")
  data = classes[name] = {'methods': [],"constructors":[],"docs": docs}



  for d in c.findall('cdecl'):
     dname = d.find('attributelist/attribute[@name="name"]').attrib["value"]
     if (d.find('attributelist/attribute[@name="kind"]').attrib["value"]!="function"): continue

     if d.find('attributelist/parmlist') is None:
       params = []
     else:
       params = [ x.attrib["value"] for x in d.findall('attributelist/parmlist/parm/attributelist/attribute[@name="type"]') ]

     rettype = d.find('attributelist/attribute[@name="type"]').attrib["value"]
     storage = getAttribute(d,"storage")

     access = getAttribute(d,"access")
     if access=="private": continue

     docs = getAttribute(d,"feature_docstring")

     data["methods"].append((dname,params,rettype,"Static" if storage=="static" else "Normal",docs))

  for d in itertools.chain(c.findall('constructor'),c.findall('extend/constructor')):
     if d.find('attributelist/parmlist') is None:
       params = []
     else:
       params = [ x.attrib["value"] for x in d.findall('attributelist/parmlist/parm/attributelist/attribute[@name="type"]') ]

     rettype = name
     access = getAttribute(d,"access")
     if access=="private": continue

     docs = getAttribute(d,"feature_docstring")
     data["methods"].append((dname,params,rettype,"Constructor",docs))

  data["bases"] = []
  for d in c.findall('attributelist/baselist/base'):
    base = d.attrib["name"]
    data["bases"].append(base)

functions = []
for d in r.findall('*//namespace/cdecl'):
  if d.find('attributelist/attribute[@name="sym_name"]') is None: continue
  if d.find('attributelist/attribute[@name="kind"]').attrib["value"]!="function": continue

  dname = d.find('attributelist/attribute[@name="sym_name"]').attrib["value"]
  if dname == "dummy": continue
  if d.find('attributelist/parmlist') is None:
    params = []
  else:
    params = [ x.attrib["value"] for x in d.findall('attributelist/parmlist/parm/attributelist/attribute[@name="type"]') ]

  rettype = getAttribute(d,"type")
  docs = getAttribute(d,"feature_docstring")

  functions.append((dname,params,rettype,docs))



treedata = {"moduleClasses": [],"moduleFunctions": [], "moduleEnums": {}}


finclude  = file(my_module+'_swiginclude.hpp','w')
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

  treedata["moduleClasses"].append({"classType": k, "classMethods": methods, "classDocs": v["docs"],"classDocslink":""})

tools = []
ioschemehelpers = []

for (name,pars,rettype,docs) in functions:
  treedata["moduleFunctions"].append({"funName": name, "funReturn": rettype, "funParams": pars, "funDocs":docs,"funDocslink":""})

enumslist = []

for k,v in enums.items():
  treedata["moduleEnums"][k] = {
    "enumDocs": v['docs'],
    "enumDocslink": "",
    "enumEntries": dict(
       (kk , {"enumEntryDocs": vv["docs"],"enumEntryDocslink":"","enumEntryVal": vv["ev"]})
          for kk,vv in v["entries"].items())
  }


treedata["moduleInheritance"] = dict((k, [k for i in v["bases"]]) for k,v in classes.items())

json.dump(treedata,file(my_module+'.json','w'),indent=True)
