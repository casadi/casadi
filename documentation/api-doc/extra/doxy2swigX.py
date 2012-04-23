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
