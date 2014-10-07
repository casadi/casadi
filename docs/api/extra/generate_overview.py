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
from casadi import *

import xml.etree.ElementTree as ET
import pydot
tree = ET.parse('XML/index.xml')
root = tree.getroot()


graph = pydot.Dot('Overview', graph_type='digraph',rankdir='LR',dpi=50,ranksep=2) 

def print_subclasses(myclass, depth=0):
  name = myclass.__name__
  refid = root.find(".//compound[name='casadi::%s']" % name).attrib["refid"]
  description = ET.parse('XML/%s.xml' % refid).getroot().findtext("doxygen/compounddef/briefdescription/para")
  graph.add_node(pydot.Node(name,URL=refid+".html",shape="box"))
  for s in myclass.__subclasses__():
    graph.add_edge(pydot.Edge(name,s.__name__,arrowtail="normal"))
    print_subclasses(s,depth=depth+1)
    
print_subclasses(Function)

graph.write_cmapx("overview.map")
graph.write_png("overview.png")
