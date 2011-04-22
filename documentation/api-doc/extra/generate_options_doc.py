"""
Grabs the pdf-compiled examples and puts them in the API tree.
"""

xml = 'XML/'
out = 'extra/'

import re
import shutil
from lxml import etree
import os
import re

# We will use an XML parser for the main index file
# The other files we will parse as text, for speed.
xmlData = etree.parse(xml+'index.xml')

classes = xmlData.findall("//compound[@kind='class']")

# The main data structure we will set up
# Keys are the full names e.g. CasADi::Sundials::CVodesIntegrator
# The values are dicts with the following structure:
# 
# parents     - the list of classes from which this class inherits
# hierarchy   - the list of ancestors of this class
# file        - the absolute path to the hpp source file
# xmlsource   - the relative path (wrt to xml) to the xml file describing the class
# hasInternal - (optional) the internal class of this class (the relative path to the xml file)
# InternalFor - the list of classes for which this class is internal
# options     - a list of dicts denoting the class's options
# monitors    - a list of dicts denoting the class's monitorables
# stats       - a list of dicts denoting the class's stats
metadata=dict()

# Text-parsing to fill in 'parents','file','hasInternal'
for c in classes:
  refid = c.attrib['refid']
  name=c.findtext("name")
  f = file(xml+refid+'.xml',"r")
  meta = dict()
  metadata[name]=meta
  meta['parents']=[]
  meta['xmlsource']=refid
  
  hasInternalAttempt=''
  fileAttemptFlag=False
  for l in f:
    if not(l.find('basecompoundref')==-1):
      m = re.search('>(.*?)<',l)
      if m:
        meta['parents'].append(m.group(1))
    if not(l.find('<location')==-1) and fileAttemptFlag:
      if not('file' in meta):
        m = re.search('file="(.*?)"',l)
        if m:
          meta['file']=m.group(1)
    if not(l.find(name)==-1) and not(l.find("<definition")==-1):
      fileAttemptFlag = True
    if not(l.find('<type><ref')==-1):
      if not('hasInternal' in meta):
        m = re.search('refid="(.*?)"',l)
        if m:
          hasInternalAttempt=m.group(1)
    if not(l.find('<name>operator-&gt;</name>')==-1):
      if not('hasInternal' in meta):
        meta['hasInternal']=hasInternalAttempt
  f.close()

# Get the parents of a class
def parents(name):
  return metadata[name]['parents'] + sum(map(parents,metadata[name]['parents']),[])

# Fill in 'hierarchy'
for k in metadata.keys():
  metadata[k]['hierarchy']=parents(k)

# Get name by specifying xml-source
def getNameByXMLsource(xmlsource):
  for name,meta in metadata.items():
    if (meta['xmlsource']==xmlsource):
      return name
  return None

# Fill in 'internalFor'
for name,meta in metadata.items():
  if ('hasInternal' in meta):
    name2 = getNameByXMLsource(meta['hasInternal'])
    if name2 is None:
      continue
    internalfor = metadata[name2]
    if ('InternalFor' in internalfor):
      internalfor['InternalFor'].append(name)
    else:
      internalfor['InternalFor']=[name]

# Inspect anything that has FXInternal as Base Class
for name,meta in metadata.items():
  if not('CasADi::FXInternal' in meta['hierarchy']):
    continue
  source = re.sub(r'\.hpp$',r'.cpp',meta['file'])
  meta['options']={}
  meta['stats']={}
  meta['monitors']={}
  f =file(source,"r")
  for l in f:
    if not(l.find('addOption')==-1):
      m = re.search(r'addOption\(\s*"(.*?)"\s*,\s*(.*?)\s*,\s*(.*?)(,\s*(.*?))?\)\s*;(\s*// (.*))',l)
      if m:
        description=[]
        if not(m.group(7) is None):
          description.append(m.group(7))
        if not(m.group(5) is None):
          description.append(m.group(5))
        meta['options'][m.group(1)]={'name': m.group(1),'type': m.group(2),'default': m.group(3),'description': '\n'.join(description), 'used': name}
    if not(l.find('stats_')==-1):
      m = re.search(r'stats_\["(.*?)"\]',l)
      if m:
        meta['stats'][m.group(1)]={'name': m.group(1), 'used': name}
    if not(l.find('monitored')==-1):
      m = re.search(r'monitored\("(.*?)"\)',l)
      if m:
        meta['monitors'][m.group(1)]={'name': m.group(1), 'used': name}
  for k in ['options','stats','monitors']:
    if len(meta[k])==0:
      del meta[k]
  f.close()
  

print metadata['CasADi::Sundials::CVodesInternal']
  
def optionsashtml(option):
  return "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>" %(option['name'],option['type'],option['default'],option['description'],option['used'])

def statsashtml(stat):
  return "<tr><td>%s</td><td>%s</td></tr>" %(stat['name'],stat['used'])

def monitorsashtml(monitor):
  return "<tr><td>%s</td><td>%s</td></tr>" %(monitor['name'],monitor['used'])

f = file(out+'options.hpp','w')

# Print out doxygen information - options
for name,meta in metadata.items():
  if not('options' in meta):
    continue

  alloptions = meta['options']
  for a in meta['hierarchy']:
    if 'options' in metadata[a]:
      alloptions.update(metadata[a]['options'])
  
  myoptionskeys = alloptions.keys()
  myoptionskeys.sort()
  
  targets = [name]
  if 'InternalFor' in meta:
    targets+=meta['InternalFor']
  for t in targets:
    f.write("/** \class %s\n" % t)
    f.write("List of available options\n")
    f.write("<table>\n")
    f.write("<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>\n")
    for k in myoptionskeys :
      f.write(optionsashtml(alloptions[k])+"\n")
    f.write( "</table>\n")
    f.write( "*/\n")

f.close()

f = file(out+'stats.hpp','w')

# Print out doxygen information - stats
for name,meta in metadata.items():
  if not('stats' in meta):
    continue

  allstats = meta['stats']
  for a in meta['hierarchy']:
    if 'stats' in metadata[a]:
      allstats.update(metadata[a]['stats'])
  
  mystatskeys = allstats.keys()
  mystatskeys.sort()
  
  targets = [name]
  if 'InternalFor' in meta:
    targets+=meta['InternalFor']
  for t in targets:
    f.write("/** \class %s\n" % t)
    f.write("List of available stats\n")
    f.write("<table>\n")
    f.write("<tr><th>Id</th><th>Used in</th></tr>\n")
    for k in mystatskeys :
      f.write(statsashtml(allstats[k])+"\n")
    f.write( "</table>\n")
    f.write( "*/\n")

f.close()

f = file(out+'monitors.hpp','w')

# Print out doxygen information - monitors
for name,meta in metadata.items():
  if not('monitors' in meta):
    continue

  allmonitors = meta['monitors']
  for a in meta['hierarchy']:
    if 'monitors' in metadata[a]:
      allmonitors.update(metadata[a]['monitors'])
  
  mymonitorskeys = allmonitors.keys()
  mymonitorskeys.sort()
  
  targets = [name]
  if 'InternalFor' in meta:
    targets+=meta['InternalFor']
  for t in targets:
    f.write("/** \class %s\n" % t)
    f.write("List of available monitors\n")
    f.write("<table>\n")
    f.write("<tr><th>Id</th><th>Used in</th></tr>\n")
    for k in mymonitorskeys :
      f.write(monitorsashtml(allmonitors[k])+"\n")
    f.write( "</table>\n")
    f.write( "*/\n")

f.close()
