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
from doxy2swig import *

import sys

import ipdb
import texttable
import re

import lxml.etree as ET
import pydot

aliases = {}
for line in file('../Doxyfile.in','r'):
  if line.startswith('ALIASES'):
    m = re.search('\+=\s*(\w+)\s*=\s*"(.*?)"',line) 
    if m:
      aliases[m.group(1)]=m.group(2)

print aliases


def astext(node,whitespace=False,escape=True,strictspacing=False):
  r = []
  if node.nodeType == node.TEXT_NODE:
    d = node.data
    if escape:
      d = d.replace('\\', r'\\\\')
      d = d.replace('"', r'\"')
    if not(whitespace):
      d = d.strip()
    r.append(d)
  elif hasattr(node,'tagName') and node.tagName=="sp":
    return " "
  elif hasattr(node,'childNodes'):
    for node in node.childNodes:
      r.append(astext(node,whitespace=whitespace,escape=escape,strictspacing=strictspacing))
  if strictspacing:
    return ("".join(r))
  else:
    return (" ".join(r)).strip()


class Doxy2SWIG_X(Doxy2SWIG):

  def __init__(self, *args, **kwargs):
    self.internal = kwargs["internal"]
    del kwargs["internal"]
    self.deprecated = kwargs["deprecated"]
    del kwargs["deprecated"]
    self.merge = kwargs["merge"]
    del kwargs["merge"]
    self.groupdoc = kwargs["groupdoc"]
    del kwargs["groupdoc"]
    Doxy2SWIG.__init__(self,*args, **kwargs)
    

    # skip sharedobjectnode stuff
    for e in self.xmldoc.getElementsByTagName("compounddef"):
      for e2 in e.getElementsByTagName("inheritancegraph"):
        for e3 in e2.getElementsByTagName("node"):
          for e4 in e3.getElementsByTagName("label"):
            if "Node" in e4.firstChild.data:
              print "skipped"
              try:
                self.xmldoc.removeChild(e)
              except:
                pass

    self.docstringmap = {}
    self.active_docstring = None
    self.add_text_counter = 0
    self.src = args[0]
    del self.ignores[self.ignores.index("programlisting")]


  def generic_parse(self, node, pad=0):
        """A Generic parser for arbitrary tags in a node.

        Parameters:

         - node:  A node in the DOM.
         - pad: `int` (default: 0)

           If 0 the node data is not padded with newlines.  If 1 it
           appends a newline after parsing the childNodes.  If 2 it
           pads before and after the nodes are processed.  Defaults to
           0.

        """
        npiece = 0
        if pad:
            npiece = self.add_text_counter
            if pad == 2:
                self.add_text('\n')                
        for n in node.childNodes:
            self.parse(n)
        if pad:
            if self.add_text_counter > npiece:
                self.add_text('\n')
                
  def clean_pieces(self, pieces):
      """Cleans the list of strings given as `pieces`.  It replaces
      multiple newlines by a maximum of 2 and returns a new list.
      It also wraps the paragraphs nicely.
      
      """
      ret = []
      count = 0
      for i in pieces:
          if i == '\n':
              count = count + 1
          else:
              if i == '";':
                  if count:
                      ret.append('\n')
              elif count > 2:
                  ret.append('\n\n')
              elif count:
                  ret.append('\n'*count)
              count = 0
              ret.append(i)

      _data = "".join(ret)
      ret = []
      for i in _data.split('\n\n'):
          if i == 'Parameters:' or i == 'Exceptions:':
              ret.extend([i, '\n-----------', '\n\n'])
          elif i.find('// File:') > -1: # leave comments alone.
              ret.extend([i, '\n'])
          else:
              if i.strip().startswith(">") or "---" in i or "===" in i or "%%newline%%" in i:
                if "%%newline%%" in i:
                  _tmp = i.replace("%%newline%%","\n")
                else:
                  _tmp = i.strip()
              else:
                _tmp = textwrap.fill(i.strip(), 80-4, break_long_words=False)
              _tmp = self.lead_spc.sub(r'\1"\2', _tmp)
              ret.extend([_tmp, '\n\n'])
      return ret
        
  def write(self, fname):
        o = my_open_write(fname)
        if self.multi:
            for p in self.pieces:
              pp = p.encode("ascii","ignore")
              pp = re.sub("[A-Z_]*_EXPORT ","",pp)
              o.write(pp)
        else:
            for p in self.clean_pieces(self.pieces):
              pp = p.encode("ascii","ignore")
              pp = re.sub("[A-Z_]*_EXPORT ","",pp)
              o.write(pp)
        o.close()
        
  def do_doxygenindex(self, node):
      self.multi = 1
      comps = node.getElementsByTagName('compound')
      for c in comps:
          refid = c.attributes['refid'].value          
          if c.attributes['kind'].value in ["example","struct","group","source","dir","file"]: continue
          filters = ["internal","interface","node"]
          filtered = False          
          for f in filters:
            if f in refid or f.capitalize() in refid: filtered = True
          if "SparsityInterface" in refid: filtered = False
          if filtered: continue

          print c.attributes['refid'].value, c.attributes['kind'].value

          fname = refid + '.xml'
          if not os.path.exists(fname):
              fname = os.path.join(self.my_dir,  fname)
          if fname.endswith("cpp.xml"): continue
          if not self.quiet:
              print "parsing file: %s"%fname
          p = Doxy2SWIG_X(fname, self.include_function_definition, self.quiet,internal=self.internal,deprecated=self.deprecated,merge=self.merge,groupdoc=self.groupdoc)
          p.generate()
          self.pieces.extend(self.clean_pieces(p.pieces))
  
  def do_hruler(self,node):
    self.add_text("\n\n" + "-"*80 + "\n\n")

  def do_heading(self,node):
    level = int(node.attributes['level'].value)-2
    text = astext(node)
    length = len(text)
    self.add_text(text+"\n")
    if level==1:
      self.add_text(("="*length) + "\n\n")
    else:
      self.add_text(("-"*length) + "\n\n")

  def do_htmlonly(self,node):
    pass
    
  def do_programlisting(self,node):
    if hasattr(node.previousSibling,'tagName') and node.previousSibling.tagName=="htmlonly" and "doctest" in node.previousSibling.firstChild.data:
      self.add_text("\n\n::\n\n")
      for codeline in node.getElementsByTagName("codeline"):
        self.add_text("  >>> " + astext(codeline,strictspacing=True).replace("\n","%%newline%%")+"\n")
    
    
  def do_verbatim(self, node):
    skipheader = False
    try:
      if node.previousSibling.tagName=="htmlonly" and "doctest" in node.previousSibling.firstChild.data:
        skipheader = True
    except:
      pass
    if not skipheader:
      self.add_text("\n\n::\n\n")
    text = node.firstChild.data
    text = "\n".join(["  " +i for i in text.split("\n")])
    
    self.add_text(text.replace("\n","%%newline%%").replace('\\', r'\\\\').replace('"', r'\"')+"\n\n")
    
  def do_table(self, node):
     
    caption = node.getElementsByTagName("caption")
    if len(caption)==1:
      self.add_text(">" + astext(caption[0]).encode("ascii","ignore")+"\n\n")
    
    rows = []
    for (i,row) in enumerate(node.getElementsByTagName("row")):
      rows.append([])
      for (j,entry) in enumerate(row.getElementsByTagName("entry")):
        rows[i].append(astext(entry,escape=False).encode("ascii","ignore").replace("&gt;",">").replace("&lt;","<"))
      
    table = texttable.Texttable(max_width=80-4)
    table.add_rows(rows)
    
    d = table.draw()
    d = d.replace('\\', r'\\\\')
    d = d.replace('"', r'\"')
    self.add_text(d)
    self.add_text("\n")
    #print table.draw()
        
    #for row in rows:
    #  self.add_text("*")
    #  r = " "
    #  for col in row:
    #    r+=col+" | "
    #  self.add_text(col[:-1]+"\n")
    #self.generic_parse(node, pad=1)

  def do_compoundname(self, node):
      self.add_text('\n\n')
      data = node.firstChild.data
      self.start_docstring(data)

  def do_compounddef(self, node):
      kind = node.attributes['kind'].value

      if kind in ('class', 'struct'):
          prot = node.attributes['prot'].value
          if prot <> 'public':
              return
          names = ('compoundname', 'briefdescription',
                   'detaileddescription', 'includes')
          first = self.get_specific_nodes(node, names)
          for n in names:
              if first.has_key(n):
                  self.parse(first[n])
          self.end_docstring()
          for n in node.childNodes:
              if n not in first.values():
                  self.parse(n)
      elif kind in ('file', 'namespace'):
          nodes = node.getElementsByTagName('sectiondef')
          for n in nodes:
              self.parse(n)
      elif kind in ('group',):
          if node.getElementsByTagName('compoundname')[0].firstChild.data.startswith("plugin_"):
            names = ('compoundname', 'briefdescription',
                     'detaileddescription', 'includes')
            first = self.get_specific_nodes(node, names)
            for n in names:
                if first.has_key(n):
                    self.parse(first[n])
            self.end_docstring()

  def do_sectiondef(self, node):
      kind = node.attributes['kind'].value
      if kind in ('public-func', 'func', 'user-defined', '','friend'):
          self.generic_parse(node)   

  def do_memberdef(self, node):
      prot = node.attributes['prot'].value
      id = node.attributes['id'].value
      kind = node.attributes['kind'].value
      tmp = node.parentNode.parentNode.parentNode
      compdef = tmp.getElementsByTagName('compounddef')[0]
      cdef_kind = compdef.attributes['kind'].value
      location = node.getElementsByTagName('location')[0].attributes['file'].value
      if prot == 'public' and not location.endswith('cpp'):
          first = self.get_specific_nodes(node, ('definition', 'name','argsstring'))
          name = first['name'].firstChild.data
          if name[:8] == 'operator': # Don't handle operators yet.
              return

          if not first.has_key('definition') or \
                 kind in ['variable', 'typedef']:
              return

          try:
            defn = first['definition'].firstChild.data + first['argsstring'].firstChild.data
          except:
            return
          target = ""
          anc = node.parentNode.parentNode
          if cdef_kind in ('file', 'namespace'):
              ns_node = anc.getElementsByTagName('innernamespace')
              if not ns_node and cdef_kind == 'namespace':
                  ns_node = anc.getElementsByTagName('compoundname')
              if ns_node:
                  ns = ns_node[0].firstChild.data
                  target = '%s::%s'%(ns, name)
              else:
                  target = '%s '%(name)
          elif cdef_kind in ('class', 'struct'):
              # Get the full function name.
              anc_node = anc.getElementsByTagName('compoundname')
              cname = anc_node[0].firstChild.data
              if kind=="friend":
                target = "friendwrap_" + name
              else:
                target = '%s::%s'%(cname, name)
              
          self.start_docstring(target,defn)
          for n in node.childNodes:
              if n not in first.values():
                  self.parse(n)
          self.end_docstring()

  def start_docstring(self,target,origin="huma kavula"):
    self.active_docstring = (target,origin)
    if target in self.docstringmap:
      self.docstringmap[target].append((origin,[]))
    else:
      self.docstringmap[target]= [(origin,[])]
    
  def end_docstring(self):
    self.active_docstring = None
  
  def add_text(self,value):
    self.add_text_counter += 1
    if self.active_docstring is None:
      self.add_text_original(value)
    else:
      #print self.active_docstring, self.docstringmap[self.active_docstring[0]], value
      if type(value) in (types.ListType, types.TupleType):
          self.docstringmap[self.active_docstring[0]][-1][1].extend(value)
      else:
          self.docstringmap[self.active_docstring[0]][-1][1].append(value)
    
  def add_text_original(self, value):
      Doxy2SWIG.add_text(self,value)
      
  def generate(self):
      try:
        publicapi = ET.parse(self.src.replace("_internal","_cluttered"))
      except Exception as e:
        publicapi = None
      Doxy2SWIG.generate(self)
      for k, v in self.docstringmap.iteritems():
        if k.startswith("plugin_"):
          groupdoc[k] = v[0][1]
          break
        # Group together
        grouped_list = []
        grouped_dict = {}
        
        def fix_signature(a): return a.replace("override","")
        
        for (origin,pieces) in v:
          origin_nostatic = origin.replace("static ","")
          m = list(re.finditer("\[DEPRECATED(:(.*?))?\]","".join(pieces)))
          deprecated = ("" if m[-1].group(2) is None else m[-1].group(2)) if len(m)>0 else None
          
          if origin == "huma kavula" and publicapi is not None:
            internal = False
          else:
            internal = True
            if publicapi is not None:
              if publicapi.xpath(".//memberdef[translate(concat(definition,argsstring),' ','')='%s' or translate(concat(definition,argsstring),' ','')='%s']" % (origin_nostatic.replace(" ",""),'static'+origin_nostatic.replace(" ",""))):
                internal = False
          if deprecated is not None: internal = False
          
          if not re.search(r"\b_\w+\(",origin_nostatic) and "(" in origin_nostatic and "~" not in origin_nostatic:
            if origin_nostatic.endswith("=0"):
              swigname = origin_nostatic[:-2]
            else:
              swigname = origin_nostatic
            # Remove return type
            
            nest=0
            sep = 0
            for i,c in enumerate(swigname):
              if c=="<" : nest+=1
              if c==">" : nest-=1
              if c == " " and nest==0: sep = i
              if c=="(" : break
            
            swigname= swigname[sep:]
          else:
            swigname = ""
            
            
          #"INTERNAL": "mark_internal(\"$decl\");"
          #"DEPRECATED": "deprecated(\"%s\");"
          #"UNSAFE": "unsafe(\"%s\");"
                            
          if internal:
            pieces = ["[INTERNAL] "] + pieces
            if "*" not in swigname and swigname not in self.internal:
            # self.internal.add("%rename(\"_internal_%s\") " + swigname + ";")
            #  print swigname
              #self.internal.add("%pythonprepend " + swigname + " %{\n _mark_internal() \n%}")
              self.internal[swigname] = "%exception " + swigname + " {\n CATCH_OR_NOT(INTERNAL_MSG() $action) \n}"
          else:
            self.internal[swigname] = ""
            
          if deprecated is not None:
            self.deprecated[swigname] = "%exception " + swigname + " {\n CATCH_OR_NOT(DEPRECATED_MSG(\"%s\") $action)\n}" % deprecated
          else:
            self.deprecated[swigname] = ""
              
          total = u"".join(pieces)
          totalnowrap = total.replace("\n"," ")
          if (aliases["noswig"] in totalnowrap) or (aliases["nopython"] in totalnowrap):
             print "skipping", origin
             continue
          if total in grouped_dict:
             grouped_dict[total][0].append(origin)
          else:
             grouped_dict[total] = ([origin],pieces)
             grouped_list.append(grouped_dict[total])
          if not self.merge:
            self.add_text_original(["%feature(\"docstring\") ", fix_signature(swigname if len(swigname)>0 else k), " \"\n\n"]+pieces+["\";\n","\n"])
            
        if self.merge:
          if len(grouped_list)==1:
            self.add_text_original(["%feature(\"docstring\") ", fix_signature(k), " \"\n"]+grouped_list[0][1]+["\";\n","\n"])
          else:
            self.add_text_original(["%feature(\"docstring\") ",fix_signature(k) , " \"\n"])
            for (origin,pieces) in grouped_list:
              if len(u"".join(pieces).rstrip())>0:
                self.add_text_original(["\n"]+["\n>  " + o.replace('"',r'\"') + '\n'  for o in origin] + ["-"*(80-8) + "\n"] + pieces + ["\n"])
            self.add_text_original(["\";\n","\n"])
  

def convert(input, output,  include_function_definition=True, quiet=False,internal=None,deprecated=None,merge=False,groupdoc=None):
    p = Doxy2SWIG_X(input, include_function_definition, quiet,internal=internal,deprecated=deprecated,merge=merge,groupdoc=groupdoc)
    p.generate()
    p.write(output)

    
if __name__ == '__main__':
    usage = __doc__
    parser = optparse.OptionParser(usage)
    parser.add_option("-n", '--no-function-definition',
                      action='store_true',
                      default=False,
                      dest='func_def',
                      help='do not include doxygen function definitions')
    parser.add_option("-q", '--quiet',
                      action='store_true',
                      default=False,
                      dest='quiet',
                      help='be quiet and minimize output')
                      
    parser.add_option("-m", '--merge',
                      action='store_true',
                      default=False,
                      dest='merge',
                      help='merge overloaded names')
    
    options, args = parser.parse_args()
    if len(args) != 4:
        parser.error("error: no input and output specified")
    
    internal = dict()
    deprecated = dict()
    groupdoc = dict()
    convert(args[0], args[1], False, options.quiet,internal=internal,deprecated=deprecated,merge=options.merge,groupdoc=groupdoc)
    file(args[2],'w').write("\n".join(sorted(filter(lambda x: len(x)>0, internal.values()))))
    file(args[3],'w').write("\n".join(sorted(filter(lambda x: len(x)>0, deprecated.values()))))
    import pickle
    filemap = pickle.load(file('filemap.pkl','r'))
    
    for k,v in groupdoc.iteritems():
      fn,n = filemap[k]
      f = file(fn.replace(".hpp","_meta.cpp"),"w")
      f.write(file('../../../misc/license_header.txt','r').read())
      f.write("""
      #include "%s"
      #include <string>

      const std::string %s::meta_doc=
      """ % (fn.split("/")[-1],n,))
      
      for p in "".join(v).split("\n"):
        f.write('"')
        f.write(p.replace('%%newline%%','\\n"\n"')+"\\n")
        f.write('"\n')
      f.write(';\n')
      f.close()
