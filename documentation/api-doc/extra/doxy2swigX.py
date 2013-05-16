#
#     This file is part of CasADi.
# 
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

def astext(node,whitespace=False,escape=True):
  r = []
  if node.nodeType == node.TEXT_NODE:
    d = node.data
    if escape:
      d = d.replace('\\', r'\\\\')
      d = d.replace('"', r'\"')
    if not(whitespace):
      d = d.strip()
    r.append(d)
  elif hasattr(node,'childNodes'):
    for node in node.childNodes:
      r.append(astext(node,whitespace=whitespace,escape=escape))
  return (" ".join(r)).strip()


class Doxy2SWIG_X(Doxy2SWIG):

  def __init__(self, *args, **kwargs):
    Doxy2SWIG.__init__(self,*args, **kwargs)
    self.docstringmap = {}
    self.active_docstring = None
    self.add_text_counter = 0


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
              if i.strip().startswith(">"):
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
              o.write(p.encode("ascii","ignore"))
        else:
            for p in self.clean_pieces(self.pieces):
              o.write(p.encode("ascii","ignore"))
        o.close()
        
  def do_doxygenindex(self, node):
      self.multi = 1
      comps = node.getElementsByTagName('compound')
      for c in comps:
          refid = c.attributes['refid'].value
          fname = refid + '.xml'
          if not os.path.exists(fname):
              fname = os.path.join(self.my_dir,  fname)
          if not self.quiet:
              print "parsing file: %s"%fname
          p = Doxy2SWIG_X(fname, self.include_function_definition, self.quiet)
          p.generate()
          self.pieces.extend(self.clean_pieces(p.pieces))
            
  def do_table(self, node):
     
    caption = node.getElementsByTagName("caption")
    if len(caption)==1:
      self.add_text(">" + astext(caption[0]).encode("ascii","ignore")+"\n")
    
    rows = []
    for (i,row) in enumerate(node.getElementsByTagName("row")):
      rows.append([])
      for (j,entry) in enumerate(row.getElementsByTagName("entry")):
        rows[i].append(astext(entry,escape=False).encode("ascii","ignore"))
      
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
      
  def do_memberdef(self, node):
      prot = node.attributes['prot'].value
      id = node.attributes['id'].value
      kind = node.attributes['kind'].value
      tmp = node.parentNode.parentNode.parentNode
      compdef = tmp.getElementsByTagName('compounddef')[0]
      cdef_kind = compdef.attributes['kind'].value
      
      if prot == 'public':
          first = self.get_specific_nodes(node, ('definition', 'name','argsstring'))
          name = first['name'].firstChild.data
          if name[:8] == 'operator': # Don't handle operators yet.
              return

          if not first.has_key('definition') or \
                 kind in ['variable', 'typedef']:
              return

          defn = first['definition'].firstChild.data + first['argsstring'].firstChild.data

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
      Doxy2SWIG.generate(self)
      for k, v in self.docstringmap.iteritems():
        # Group together
        grouped_list = []
        grouped_dict = {}
        for (origin,pieces) in v:
          total = u"".join(pieces)
          if total in grouped_dict:
             grouped_dict[total][0].append(origin)
          else:
             grouped_dict[total] = ([origin],pieces)
             grouped_list.append(grouped_dict[total])
        if len(grouped_list)==1:
          self.add_text_original(["%feature(\"docstring\") ", k, " \"\n"]+grouped_list[0][1]+["\";\n","\n"])
        else:
          self.add_text_original(["%feature(\"docstring\") ",k , " \"\n"])
          for (origin,pieces) in grouped_list:
            if len(u"".join(pieces).rstrip())>0:
              self.add_text_original(["\n"]+["\n>  " + o + '\n'  for o in origin] + ["-"*(80-8) + "\n"] + pieces + ["\n"])
          self.add_text_original(["\";\n","\n"])
      
def convert(input, output, include_function_definition=True, quiet=False):
    p = Doxy2SWIG_X(input, include_function_definition, quiet)
    p.generate()
    p.write(output)

def main():
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
    
    options, args = parser.parse_args()
    if len(args) != 2:
        parser.error("error: no input and output specified")

    convert(args[0], args[1], False, options.quiet)
    

if __name__ == '__main__':
    main()
