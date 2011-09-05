from casadi import SX, SXMatrix, getOperatorRepresentation

import pydot
      
def adjSXDot(graph,adj):
  """
  Given an adjacency dict structure and a graph, add nodes an edges to the graph to represent an SX tree
  """
  for k,v in adj.iteritems():
    if k.isLeaf():
      style = "solid"
      if k.isConstant():
        style = "bold"
      graph.add_node(pydot.Node(str(k.__hash__()),label=str(k),shape="box",style=style))
    else:
      if not(k.isCommutative()):
        s = getOperatorRepresentation(k,["| <f0> | ", " | <f1> |"])
        if s.startswith("(|") and s.endswith("|)"):
          s=s[2:-2]
        
        graph.add_node(pydot.Node(str(k.__hash__()),label=s,shape='Mrecord'))
        for i,n in enumerate(v):
          graph.add_edge(pydot.Edge(str(n.__hash__()),str(k.__hash__())+":f%d" % i))
      else:
        s = getOperatorRepresentation(k,[".", "."])
        if s.startswith("(.") and s.endswith(".)"):
          s=s[2:-2]
        if s.startswith("(") and s.endswith(")"):
          s=s[1:-1]
        graph.add_node(pydot.Node(str(k.__hash__()),label=s,shape='oval'))
        for i,n in enumerate(v):
          graph.add_edge(pydot.Edge(str(n.__hash__()),str(k.__hash__())))
   
def populateAdj(s,adj):
  """ This function adds to the adjacency dict recursively.
      Note: This is certainly not efficient, and will result in stack overflow for large graphs.
            However, there's not much point anyway in visualising large graphs
  """
  adj[s] = []
  if not(s.isLeaf()):
    for i in range(s.getNdeps()):
      populateAdj(s.getDep(i),adj)
      adj[s].append(s.getDep(i))
      
      
def invAdj(adj):
  """
  Returns an inverted version of the adjacency dict mapping
  """
  inv = {}
  for k,v in adj.iteritems():
    for n in v:
      if not(n in inv):
        inv[n] = [k]
      else:
        inv[n].append(k)
  return inv
      
def dotgraph(s,direction="BT"):
  """
  direction   one of "BT", "LR", "TB", "RL"
  """
      
  graph = pydot.Dot('G', graph_type='digraph',rankdir=direction)
    
  if isinstance(s,SXMatrix):
    sp = s.sparsity()
    adj = {}
    row = sp.getRow()
    for k in range(s.size()):
      sx = s.at(k)
      i = row[k]
      j = sp.col()[k]
      graph.add_edge(pydot.Edge(str(sx.__hash__()),"Matrix:f%d%d" % (i,j)))
      populateAdj(sx,adj)
    invadj = invAdj(adj)
    adjSXDot(graph,adj)
    
    label = '<<TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0">'
    for i in range(s.size1()):
      label+="<TR>"
      for j in range(s.size2()):
        k = sp.getNZ_const(i,j)
        if k==-1:
          label+="<TD>.</TD>"
        else:
          if not(s.at(k) in invadj) and s.at(k).isLeaf():
            graph.del_node(str(s.at(k).__hash__()))
            graph.del_edge(str(s.at(k).__hash__()),dst="Matrix:f%d%d" % (i,j))
            label+="<TD>%s</TD>" % str(s.at(k))
          else:
            label+="<TD PORT='f%d%d'> <I>(%d,%d)</I> </TD>" % (i,j,i,j)
      label+="</TR>"
    label+="</TABLE>>"
    graph.add_node(pydot.Node("Matrix",label=label,shape='plaintext'))
    
  elif isinstance(s,SX):
    adj = {}
    populateAdj(s,adj)
    adjSXDot(graph,adj)
  
  return graph


def dotdraw(s,direction="RL"):
  """
  direction   one of "BT", "LR", "TB", "RL"
  """

  from pylab import imread, imshow,show,figure, axes
  if hasattr(show,'__class__') and show.__class__.__name__=='PylabShow':  # catch pyreport
    figure_name = '%s%d.%s' % ( show.basename, len(show.figure_list), show.figure_extension )
    show.figure_list += (figure_name, )
    dotgraph(s,direction=direction).write_pdf(figure_name)
    print "Here goes figure %s (dotdraw)" % figure_name
  else:
    temp="_temp.png"
    dotgraph(s,direction=direction).write_png(temp)
    im = imread(temp)
    figure()
    ax = axes([0,0,1,1], frameon=False)
    ax.set_axis_off()
    imshow(im)
    show()
