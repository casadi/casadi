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
