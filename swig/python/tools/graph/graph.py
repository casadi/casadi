from casadi import SX, SXMatrix, getOperatorRepresentation

try:
  import pydot
except:
  raise Exception("To use the functionality of casadi.tools.graph, you need to have pydot Installed. Try `easy_install pydot`.")  
      
def adjSXDot(graph,adj):
  """
  Given an adjacency dict structure and a pydot graph, adds nodes an edges to the graph to represent an SX tree
  """
  for k,v in adj.iteritems():
    if k.isLeaf():
      style = "solid" # Symbolic nodes are represented box'es
      if k.isConstant():
        style = "bold" # Constants are represented by bold box'es
      graph.add_node(pydot.Node(str(k.__hash__()),label=str(k),shape="box",style=style))
    else:
      if not(k.isCommutative()):
        # Non-commutative operators are represented by 'record' shapes.
        # The dependencies have different 'ports' where arrows should arrive.
        s = getOperatorRepresentation(k,["| <f0> | ", " | <f1> |"])
        if s.startswith("(|") and s.endswith("|)"):
          s=s[2:-2]
        
        graph.add_node(pydot.Node(str(k.__hash__()),label=s,shape='Mrecord'))
        for i,n in enumerate(v):
          graph.add_edge(pydot.Edge(str(n.__hash__()),str(k.__hash__())+":f%d" % i))
      else: 
       # Commutative operators can be represented more compactly as 'oval' shapes.
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
  Creates and returns a pydot graph structure that represents an SX or SXMatrix.
  
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
    
    # The Matrix grid is represented by a html table with 'ports'
    label = '<<TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0">'
    for i in range(s.size1()):
      label+="<TR>"
      for j in range(s.size2()):
        k = sp.getNZ_const(i,j)
        if k==-1:
          label+="<TD>.</TD>"
        else:
          if not(s.at(k) in invadj) and s.at(k).isLeaf():
            # Compaction rule: if a leaf node points to matrix entry port,
            # and does not point anywhere else, then we remove node and edge
            # and display the leaf's reprentation directly into the matrix.
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


def dotsave(s,format='ps',filename="temp",direction="RL"):
  """
  Make a drawing of an SX or SXMatrix and save it.
  
  format can be one of:
    dot canon cmap cmapx cmapx_np dia dot fig gd gd2 gif hpgl imap imap_np
    ismap jpe jpeg jpg mif mp pcl pdf pic plain plain-ext png ps ps2 raw
    svg svgz vml vmlz vrml vtx wbmp xdot xlib
    
  direction   one of "BT", "LR", "TB", "RL"
  
  """
  g = dotgraph(s,direction=direction)
  if hasattr(g,'write_'+format):
    getattr(g,'write_'+format)(filename)
  else:
    s = "Unknown format '%s'. Please pick one of the following:\n" % format
    l = ['dot']
    for n in dir(g):
      if n.startswith("write_"):
        l.append(n[6:])
    s+= " ".join(l)
    raise Exception(s)
  
def dotdraw(s,direction="RL"):
  """
  Make a drawing of an SX or SXMatrix and display it.
  
  direction   one of "BT", "LR", "TB", "RL"
  """
  
  try:  # Check if we have pylab
    from pylab import imread, imshow,show,figure, axes
  except:
    # We don't have pylab, so just write out to file
    print "casadi.tools.graph.dotdraw: no pylab detected, will not show drawing on screen."
    dotgraph(s,direction=direction).write_ps("temp.ps")
    return
   
  if hasattr(show,'__class__') and show.__class__.__name__=='PylabShow':  
    # catch pyreport case, so we have true vector graphics
    figure_name = '%s%d.%s' % ( show.basename, len(show.figure_list), show.figure_extension )
    show.figure_list += (figure_name, )
    dotgraph(s,direction=direction).write_pdf(figure_name)
    print "Here goes figure %s (dotdraw)" % figure_name
  else:
    # Matplotlib does not allow to display vector graphics on screen, 
    # so we fall back to png
    temp="_temp.png"
    dotgraph(s,direction=direction).write_png(temp)
    im = imread(temp)
    figure()
    ax = axes([0,0,1,1], frameon=False)
    ax.set_axis_off()
    imshow(im)
    show()
