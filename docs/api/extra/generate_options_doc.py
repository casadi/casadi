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
"""
This script reads XML doxygen output and reads source-files to harvest Options.

This script generates a header file with option information that doxygen can use.
"""


xml = 'XML_internal/'
out = 'extra/'

import re
import shutil
from lxml import etree
import os
import re
import fnmatch

from generate_options_helpers import addExtra

# We will use an XML parser for the main index XML file
xmlData = etree.parse(xml+'index.xml')

classes = xmlData.findall("//compound[@kind='class']")

# The main data structure we will set up
# Keys are the full names e.g. casadi::Sundials::CVodesIntegrator
# The values are dicts with the following structure:
# 
# parents     - the list of classes from which this class inherits
# hierarchy   - the list of ancestors of this class. First element is the parent.
# file        - the absolute path to the hpp source file
# xmlsource   - the relative path (wrt to xml) to the xml file describing the class
# hasInternal - (optional) the internal class of this class (the relative path to the xml file)
# InternalFor - the list of classes for which this class is internal
# options     - a list of dicts denoting the class's options
# monitors    - a list of dicts denoting the class's monitorables
# stats       - a list of dicts denoting the class's stats
# inputscheme  - the input scheme of the function
# outputscheme - the output scheme of the function
metadata=dict()

# Text-parsing to fill in 'parents','file','hasInternal'
for c in classes:
  refid = c.attrib['refid']
  name=c.findtext("name")
  
  constructor_name = name + "::" + name.split("::")[-1]
  
  # Parse the XML file that describes the class
  f = etree.parse(xml+refid+'.xml')
  
  meta = dict()
  metadata[name]=meta
  meta['parents']=[]
  meta['xmlsource']=refid
  
  if f.getroot().find("compounddef/briefdescription/para") is not None:
    meta["brief"] = "".join(list(f.getroot().find("compounddef/briefdescription/para").itertext())).strip()
  else:
    meta["brief"] = ""
 
  # find parents
  for e in f.getroot().findall("compounddef/basecompoundref"):
    if e.text.find("std::")==-1: # Exclude standard library from hierarchy
      meta['parents'].append(e.text)

  # find the internal class of this class, if any
  temp = f.find("//memberdef[name='operator->']/type/ref")
  if not(temp is None):
    meta['hasInternal']=temp.attrib["refid"]
  
  # find the file location of this class
  temp = f.find("//memberdef[definition='%s']/location" % constructor_name)
  if not(temp is None):
    meta['file']=temp.attrib["file"]
    
print metadata["casadi::Integrator"]

for v in metadata.values(): # Remove templating if not in index.xml
  for i,vv in enumerate(v['parents']):
    if not(vv in metadata):
      v['parents'][i] = vv.split("<")[0]
      
# Get the parents of a class
def parents(name):
  if name not in metadata:
    return []
  return metadata[name]['parents'] + sum(map(parents,metadata[name]['parents']),[])

# Fill in 'hierarchy'
for k in metadata.keys():
  metadata[k]['hierarchy']=parents(k)

# Get name by specifying xml-source
def nameByXMLsource(xmlsource):
  for name,meta in metadata.items():
    if (meta['xmlsource']==xmlsource):
      return name
  return None

# Fill in 'internalFor'
for name,meta in metadata.items():
  if ('hasInternal' in meta):
    name2 = nameByXMLsource(meta['hasInternal'])
    if name2 is None:
      continue
    internalfor = metadata[name2]
    if ('InternalFor' in internalfor):
      internalfor['InternalFor'].append(name)
    else:
      internalfor['InternalFor']=[name]

# http://code.activestate.com/recipes/499305/
def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below
    supplied root directory.'''
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)
  
# locate the namespacecasadi.xml
xmlNS = etree.parse(list(locate('namespacecasadi.xml',root=xml))[0])

# construct a table with enum info
enums = {}

for enum in xmlNS.findall("//memberdef[@kind='enum']"):
  name=enum.findtext("name")
  enums[name]=[]
  description = ""
  
  if enum.findtext("briefdescription/para"):
    description+="".join(list(enum.find("briefdescription/para").itertext())).strip()
  if enum.findtext("detaileddescription/para"):
    description+= "\n" + "".join(list(enum.find("detaileddescription/para").itertext())).strip()
  enums[name].append(description)
  
  for enumvalue in enum.findall("enumvalue"):
    name_=enumvalue.findtext("name")
    description = ""
    if enumvalue.findtext("briefdescription/para"):
      description += enumvalue.findtext("briefdescription/para").strip()
    if enumvalue.findtext("detaileddescription/para"):
      description += "\n" + enumvalue.findtext("detaileddescription/para").strip()
    enums[name].append({'name': name_, 'description': description.strip()})
    
from pyparsing import *


multiline_string = OneOrMore(dblQuotedString)

def removequotes(s):
  return s[1:-1]
  
parse_quoted_string = multiline_string.copy().setParseAction(lambda x: "".join(map(removequotes,x)))
parse_quoted_string_keep = multiline_string.copy().setParseAction(lambda x: '"' + ("".join(map(removequotes,x)))+'"')
parse_type = Word( srange("[A-Z]") + "_").setResultsName("type")

comma = Suppress(Literal(","))

def brackets(a):
  return Suppress(Literal("{")) + Optional(a) + Suppress(Literal("}"))

def gbrackets(a):
  return Group(brackets(a))
  
def delimitedListTrailing( expr, delim=","):
    return ( expr + ZeroOrMore( Suppress( delim ) + expr ) + Optional( Suppress(delim) ) )


parse_optiondata = parse_type + comma + parse_quoted_string
parse_optiondata_collection = gbrackets(parse_optiondata)

parse_option = parse_quoted_string + comma + parse_optiondata_collection
parse_option_collection = gbrackets(parse_option)

parse_options = delimitedListTrailing(parse_option_collection, delim=",")
parse_options_collection = gbrackets(parse_options)

parse_options_group = (brackets( Suppress(Literal("&")) + Word(alphanums + ".:-_"))) | parse_options_collection
parse_options_group_collection = Optional(Suppress(Literal("="))) + brackets(delimitedList(parse_options_group)) + Suppress(Literal(";"))
 
# Inspect anything that has FXInternal as Base Class
for name,meta in metadata.items():
  if not('file' in meta):
    meta['options']={}
    meta['stats']={}
    meta['monitors']={}
    meta['optionproviders'] = []
    continue
  source = re.sub(r'(_impl)?\.hpp$',r'.cpp',meta['file'])
  meta['options']={}
  meta['stats']={}
  meta['monitors']={}
  meta['optionproviders'] = []
  try:
    f =file(source,"r")
  except:
    continue
  lines = f.readlines()
  linec = 0
  while linec<len(lines):
    l = lines[linec]
    linec+=1
    if 'Options' in l:
      if re.search('Options .*::options_',l):
        start = linec
        while "};" not in l:
          l = lines[linec]
          linec+=1
        stop = linec
        
        optionsdict = "".join(lines[start:stop])

        results = parse_options_group_collection.parseString(optionsdict)

        for optiongroup in results:
          if isinstance(optiongroup,str):
            m = re.search("(.*?)::options_",optiongroup)
            assert m
            op = m.group(1)
            if "casadi" not in op: op = "casadi::" + op
            meta['optionproviders'].append(op)
            continue
          if len(optiongroup[0])==0:
            continue
          for op_name, op_data in optiongroup:   
            d = meta['options'][op_name]={'name': op_name,"type": op_data[0],'used': name,'default':'','description':op_data[1],'inherit': False}
      
    if not(l.find('ops_')==-1):
      m = re.search(r'ops_\["(.*?)"\]\s*=\s*(.*?);',l)
      if m:
        meta['options'][m.group(1)]={'name': m.group(1),'type': m.group(2),'default': '','description': '', 'used': name}
    if not(l.find('stats_')==-1):
      m = re.search(r'stats_\["(.*?)"\]',l)
      if m:
        meta['stats'][m.group(1)]={'name': m.group(1), 'used': name}
    if not(l.find('INPUTSCHEME')==-1):
      m = re.search(r'INPUTSCHEME\((.*?)\)',l)
      if m:
        meta['inputscheme']=m.group(1)
    if not(l.find('OUTPUTSCHEME')==-1):
      m = re.search(r'OUTPUTSCHEME\((.*?)\)',l)
      if m:
        meta['outputscheme']=m.group(1)

  for k in ['options','stats','monitors']:
    if len(meta[k])==0:
      del meta[k]
  f.close()

addExtra(metadata)

def htmlescape(h):
  return h.replace(">","&gt;").replace("<","&lt;")

def newline2br(a):
  return a.replace("\n","<br />")
  
def optionsashtml(option,used=True):
  if used:
    return "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>" %(option['name'],option['type'],newline2br(htmlescape(option['description'])),newline2br(option['used']))
  else:
    return "<tr><td>%s</td><td>%s</td><td>%s</td></tr>" %(option['name'],option['type'],newline2br(htmlescape(option['description'])))

def statsashtml(stat,used=True):
  if used:
    return "<tr><td>%s</td><td>%s</td></tr>" %(stat['name'],stat['used'])
  else:
    return "<tr><td>%s</td></tr>" % stat['name']

def monitorsashtml(monitor,used=True):
  if used:
    return "<tr><td>%s</td><td>%s</td></tr>" %(monitor['name'],monitor['used'])
  else:
    return "<tr><td>%s</td></tr>" % monitor['name']


def update_no_overwrite(orig,new):
  for (k,v) in new.iteritems():
    if not(k in orig):
      orig[k] = v
        
def update_overwrite(orig,new):
  import copy
  for (k,v) in new.iteritems():
    if not(k in orig):
      orig[k] = copy.copy(v)
    else:
      if "inherit" in v and v["inherit"]:
        orig[k]["description"] = orig[k]["description"] + "\n" +  v["description"]
        orig[k]["used"] = orig[k]["used"] + "\n" + v["used"]
      else:
        orig[k] = copy.copy(v)
        
f = file(out+'b0_options.hpp','w')

#raise Exception("")
fdiagram = file(out+'e0_diagram.hpp','w')

filemap = {}
for name,meta in sorted(metadata.items()):
  if "casadi::PluginInterface" in meta["hierarchy"] and 'casadi::FunctionInternal' not in meta["parents"]: 
    m = re.search("'(\w+)' plugin for (\w+)",meta["brief"])
    if m:
      filemap["plugin_%s_%s" % ( m.group(2),m.group(1))] = (meta["file"],name)

import pickle
pickle.dump(filemap,file(out+"filemap.pkl","w"))

# Print out doxygen information - options
for name,meta in sorted(metadata.items()):
  if not('options' in meta):
    meta['options'] = {}

  optionproviders = [meta['options']]
  for a in meta['optionproviders']:
    if a in metadata and 'options' in metadata[a]:
      optionproviders.append(metadata[a]['options'])
  
  alloptions = {}
  
  for optionprovider in reversed(optionproviders):
      update_overwrite(alloptions,optionprovider)
      
  myoptionskeys = alloptions.keys()
  if len(myoptionskeys)==0:
    continue
  myoptionskeys.sort()
  
  for t in [name]:
    f.write("/// \cond INTERNAL\n")
    f.write("/** \class %s\n\\n\n\\par\n" % t)
    f.write("<a name='options'></a><table>\n")
    f.write("<caption>List of available options</caption>\n")
    f.write("<tr><th>Id</th><th>Type</th><th>Description</th><th>Used in</th></tr>\n")
    for k in myoptionskeys :
      f.write(optionsashtml(alloptions[k])+"\n")
    f.write( "</table>\n")
    f.write( "*/\n")
    f.write("/// \endcond\n")
    
    fdiagram.write("/// \cond INTERNAL\n")
    fdiagram.write("/** \class %s\n" % t)
    fdiagram.write("""<a name="diagram" id="diagram"></a><h2>Diagrams</h2>""")
    fdiagram.write( "*/\n")
    fdiagram.write("/// \endcond\n")
    

  if "casadi::PluginInterface" in meta["hierarchy"]:
    m = re.search("'(\w+)' plugin for (\w+)",meta["brief"])
    if not m:
      print "This plugin is undocumented. add \\pluginbrief{class,name} to it: " + meta["file"]
      #print meta["file"]
    else:
      f.write("/** \\addtogroup plugin_%s_%s\n\\n\n\\par\n" % (m.group(2),m.group(1)))
      f.write("<a name='options'></a><table>\n")
      f.write("<caption>List of available options</caption>\n")
      f.write("<tr><th>Id</th><th>Type</th><th>Description</th></tr>\n")
      for k in myoptionskeys :
        try:
          #if alloptions[k]["used"] in metadata and "casadi::PluginInterface" in metadata[alloptions[k]["used"]]["parents"]:
          if alloptions[k]["used"] in metadata and "casadi::PluginInterface" in metadata[metadata[alloptions[k]["used"]]["parents"][0]]["hierarchy"]:
            f.write(optionsashtml(alloptions[k],used=False)+"\n")
        except:
          pass
      f.write( "</table>\n")
      f.write( "*/\n")
    
    t = name
    print "test", name, myoptionskeys
    f.write("/** \\addtogroup general_%s\n\\n\n\\par\n" % t.replace("casadi::","") )
    f.write("<a name='options'></a><table>\n")
    f.write("<caption>List of available options</caption>\n")
    f.write("<tr><th>Id</th><th>Type</th><th>Description</th><th>Used in</th></tr>\n")
    for k in myoptionskeys :
      f.write(optionsashtml(alloptions[k])+"\n")
    f.write( "</table>\n")
    f.write( "*/\n")
    
    fdiagram.write("/** \class %s\n" % t)
    fdiagram.write("""<a name="diagram" id="diagram"></a><h2>Diagrams</h2>""")
    fdiagram.write( "*/\n")
  if "FunctionInternal" in name:
    f.write("/** \class casadi::Function\n")
    f.write("<a name='options'></a><table>\n")
    f.write("<caption>List of available options</caption>\n")
    f.write("<tr><th>Id</th><th>Type</th><th>Description</th><th>Used in</th></tr>\n")
    for k in myoptionskeys :
      f.write(optionsashtml(alloptions[k])+"\n")
    f.write( "</table>\n")
    f.write( "*/\n")

f.close()

"""
f = file(out+'d0_stats.hpp','w')

# Print out doxygen information - stats
for name,meta in sorted(metadata.items()):
  if not('stats' in meta):
    continue

  allstats = meta['stats']
  for a in meta['hierarchy']:
    if 'stats' in metadata[a]:
      update_no_overwrite(allstats,metadata[a]['stats'])
  
  mystatskeys = allstats.keys()
  mystatskeys.sort()
  if len(mystatskeys)==0:
    continue
    
  for t in [name]:
    f.write("/// \cond INTERNAL\n")
    f.write("/** \class %s\n\\n\n\\par\n" % t)
    f.write("<a name='stats'></a><table>\n")
    f.write("<caption>List of available stats</caption>\n")
    f.write("<tr><th>Id</th><th>Used in</th></tr>\n")
    for k in mystatskeys :
      f.write(statsashtml(allstats[k])+"\n")
    f.write( "</table>\n")
    f.write( "*/\n")
    f.write("/// \endcond\n")

  if "casadi::PluginInterface" in meta["hierarchy"] and 'casadi::FunctionInternal' not in meta["parents"]:
        
    m = re.search("'(\w+)' plugin for (\w+)",meta["brief"])
    if not m:
      print "This plugin is undocumented. add \\pluginbrief{class,name} to it: " + meta["file"]
      #print meta["file"]
    else:
      f.write("/** \\addtogroup plugin_%s_%s\n\\n\n\\par\n" % (m.group(2),m.group(1)))
      f.write("<a name='stats'></a><table>\n")
      f.write("<caption>List of available stats</caption>\n")
      f.write("<tr><th>Id</th></tr>\n")
      for k in mystatskeys :
        if allstats[k]["used"] == name:
          f.write(statsashtml(allstats[k],used=False)+"\n")
      f.write( "</table>\n")
      f.write( "*/\n")
      
  if 'InternalFor' in meta:
    for t in meta['InternalFor']:
      if "casadi::PluginInterface" in meta["parents"]:
        f.write("/** \\addtogroup general_%s\n\\n\n\\par\n" % t.replace("casadi::","") )
      else:
        f.write("/** \class %s\n\\n\n\\par\n" % t)
      f.write("<a name='stats'></a><table>\n")
      f.write("<caption>List of available stats</caption>\n")
      f.write("<tr><th>Id</th><th>Used in</th></tr>\n")
      for k in mystatskeys :
        f.write(statsashtml(allstats[k])+"\n")
      f.write( "</table>\n")
      f.write( "*/\n")
    
f.close()
"""

"""
f = file(out+'c0_monitors.hpp','w')

# Print out doxygen information - monitors
for name,meta in sorted(metadata.items()):
  if not('monitors' in meta):
    continue

  allmonitors = meta['monitors']
  for a in meta['hierarchy']:
    if 'monitors' in metadata[a]:
      update_no_overwrite(allmonitors,metadata[a]['monitors'])
  
  mymonitorskeys = allmonitors.keys()
  mymonitorskeys.sort()
  
  if len(mymonitorskeys)==0:
    continue
    
  for t in [name]:
    f.write("/// \cond INTERNAL\n")
    f.write("/** \class %s\n\\n\n\\par\n" % t)
    f.write("<a name='monitors'></a><table>\n")
    f.write("<caption>List of available monitors</caption>\n")
    f.write("<tr><th>Id</th><th>Used in</th></tr>\n")
    for k in mymonitorskeys :
      f.write(monitorsashtml(allmonitors[k])+"\n")
    f.write( "</table>\n")
    f.write( "*/\n")
    f.write("/// \endcond\n")

  if "casadi::PluginInterface" in meta["hierarchy"] and 'casadi::FunctionInternal' not in meta["parents"]:
        
    m = re.search("'(\w+)' plugin for (\w+)",meta["brief"])
    if not m:
      print "This plugin is undocumented. add \\pluginbrief{class,name} to it: " + meta["file"]
      #print meta["file"]
    else:
      f.write("/** \\addtogroup plugin_%s_%s\n\\n\n\\par\n" % (m.group(2),m.group(1)))
      f.write("<a name='monitors'></a><table>\n")
      f.write("<caption>List of available monitors</caption>\n")
      f.write("<tr><th>Id</th></tr>\n")
      for k in mymonitorskeys :
        if allmonitors[k]["used"] == name:
          f.write(monitorsashtml(allmonitors[k],used=False)+"\n")
      f.write( "</table>\n")
      f.write( "*/\n")
      
  if 'InternalFor' in meta: 
    for t in meta['InternalFor']:
      if "casadi::PluginInterface" in meta["parents"]:
        f.write("/** \\addtogroup general_%s\n\\n\n\\par\n" % t.replace("casadi::","") )
      else:
        f.write("/** \class %s\n\\n\n\\par\n" % t)
      f.write("<a name='monitors'></a><table>\n")
      f.write("<caption>List of available monitors</caption>\n")
      f.write("<tr><th>Id</th><th>Used in</th></tr>\n")
      for k in mymonitorskeys :
        f.write(monitorsashtml(allmonitors[k])+"\n")
      f.write( "</table>\n")
      f.write( "*/\n")

f.close()
"""


f = file(out+'a0_schemes.hpp','w')

def enumsashtml(n,title):
  s=""
  if (n in enums):
    s+= "<a name='schemes'></a><table>\n"
    
    num = ""
    
    for i in range(len(enums[n])-1):
      m = enums[n][1+i]
      if re.search(r'_NUM_IN$',m['name']):
        num = m['name']
      if re.search(r'_NUM_OUT$',m['name']):
        num = m['name']
    s+= "<caption>%s  (%s = %d)</caption>\n" % (title,num,len(enums[n])-2)
    s+= "<tr><th>Full name</th><th>Short</th><th>Description</th></tr>\n"
    for i in range(len(enums[n])-1):
      m = enums[n][1+i]
      if re.search(r'_NUM_IN$',m['name']):
        continue
      if re.search(r'_NUM_OUT$',m['name']):
        continue
      if m['description']=='': continue
      shorthand = "_".join(m['name'].split("_")[1:]).lower()
      s+="<tr><td>%s</td><td>%s</td><td>%s</td></tr>\n" % (m['name'],shorthand,m['description'])
    s+="</table>\n"
  return s

# Write out all input/output information
for scheme in enums:
  types = ["input","output","struct"]
  for t in types:
    if re.search(t,scheme,re.IGNORECASE):
      f.write("/** \defgroup scheme_%s\n" % scheme)
      f.write(enumsashtml(scheme,"%s scheme: casadi::%s" % (t.capitalize(),scheme)))
      f.write( "*/\n")

# Print out doxygen information - schemes
for name,meta in sorted(metadata.items()):

  inputscheme = None
  outputscheme = None
  if not('inputscheme' in meta):
    for a in meta['hierarchy']:
      if a in metadata and 'inputscheme' in metadata[a]:
        inputscheme = metadata[a]['inputscheme']
  else:
    inputscheme = meta['inputscheme']

  if not('outputscheme' in meta):
    for a in meta['hierarchy']:
      if  a in metadata and 'outputscheme' in metadata[a]:
        outputscheme = metadata[a]['outputscheme']
  else:
    outputscheme = meta['outputscheme']
  
  targets = [name]
  if 'InternalFor' in meta:
    targets+=meta['InternalFor']
    
  for t in [name]:
    if not(inputscheme is None) or not(outputscheme is None):
      f.write("/// \cond INTERNAL\n")
      f.write("/** \class %s\n\\n\n\\par\n" % t)
      if not(inputscheme is None) and not(outputscheme is None):
        f.write("@copydoc scheme_%s\n" % inputscheme)
        f.write("<br/>\n")
        f.write("@copydoc scheme_%s\n" % outputscheme)
      elif outputscheme is None:
        f.write("@copydoc scheme_%s\n" % inputscheme)
      elif inputscheme is None:
        f.write("@copydoc scheme_%s\n" % outputscheme)
      f.write( "*/\n")
      f.write("/// \endcond\n")
      
  if 'InternalFor' in meta:
    for t in meta['InternalFor']:
      if not(inputscheme is None) or not(outputscheme is None):
        if "casadi::PluginInterface" in meta["parents"]:
          f.write("/** \\addtogroup general_%s\n\\n\n\\par\n" % t.replace("casadi::","") )
        else:
          f.write("/** \class %s\n\\n\n\\par\n" % t)
        if not(inputscheme is None) and not(outputscheme is None):
          f.write("@copydoc scheme_%s\n" % inputscheme)
          f.write("<br/>\n")
          f.write("@copydoc scheme_%s\n" % outputscheme)
        elif outputscheme is None:
          f.write("@copydoc scheme_%s\n" % inputscheme)
        elif inputscheme is None:
          f.write("@copydoc scheme_%s\n" % outputscheme)
        f.write( "*/\n")
   
f.close()
