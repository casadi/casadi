from casadi import SX, SXMatrix, MX, getOperatorRepresentation

try:
  import pydot
except:
  raise Exception("To use the functionality of casadi.tools.graph, you need to have pydot Installed. Try `easy_install pydot`.")  
      
#import ipdb
 
def addDependency(master,slave,dep={},invdep={}):
  #print master.__hash__() , " => ", slave.__hash__(), "   ", master , " => ", slave
  if master in dep:
    dep[master].add(slave)
  else:
    dep[master]= set([slave])
  if slave in invdep:
    invdep[slave].add(master)
  else:
    invdep[slave]= set([master])
    
def addDependencies(master,slaves,dep={},invdep={}):
  for slave in slaves:
    addDependency(master,slave,dep = dep,invdep = invdep)
  for slave in slaves:
    dependencyGraph(slave,dep = dep,invdep = invdep)

def dependencyGraph(s,dep = {},invdep = {}):
  if isinstance(s,SXMatrix):
    addDependencies(s,list(s.data()),dep = dep,invdep = invdep)
  elif isinstance(s,SX):
    if not(s.isLeaf()):
      addDependencies(s,[s.getDep(k) for k in range(s.getNdeps())],dep = dep,invdep = invdep)
  elif isinstance(s,MX):
    addDependencies(s,[s.getDep(k) for k in range(s.getNdeps())],dep = dep,invdep = invdep)
  return (dep,invdep)
  
class DotArtist:
  def __init__(self,s,dep={},invdep={},graph=None,artists={}):
    self.s = s
    self.dep = dep
    self.invdep = invdep
    self.graph = graph
    self.artists = artists

class MXSymbolicArtist(DotArtist):
  def draw(self):
    s = self.s
    graph = self.graph
    sp = s.sparsity()
    row = sp.getRow()
      
    # The Matrix grid is represented by a html table with 'ports'
    label = '<<TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0">'
    label+="<TR><TD COLSPAN='%d'>%s: <font color='#666666'>%s</font></TD></TR>" % (s.size2(),s.getName(), s.dimString())
    for i in range(s.size1()):
      label+="<TR>"
      for j in range(s.size2()):
        k = sp.getNZ_const(i,j)
        if k==-1:
          label+="<TD>.</TD>"
        else:
          label+="<TD PORT='f%d%d' BGCOLOR='#eeeeee'> <font color='#666666'>(%d,%d)</font> </TD>" % (i,j,i,j)
      label+="</TR>"
    label+="</TABLE>>"
    graph.add_node(pydot.Node(str(self.s.__hash__()),label=label,shape='plaintext'))

class MXConstantArtist(DotArtist):
  def draw(self):
    s = self.s
    graph = self.graph
    sp = s.sparsity()
    row = sp.getRow()
    M = s.getConstant()
    # The Matrix grid is represented by a html table with 'ports'
    label = '<<TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0">'
    label+="<TR><TD COLSPAN='%d'><font color='#666666'>%s</font></TD></TR>" % (s.size2(), s.dimString())
    for i in range(s.size1()):
      label+="<TR>"
      for j in range(s.size2()):
        k = sp.getNZ_const(i,j)
        if k==-1:
          label+="<TD>.</TD>"
        else:
          label+="<TD PORT='f%d%d' BGCOLOR='#eeeeee'> %s </TD>" % (i,j,M[i,j])
      label+="</TR>"
    label+="</TABLE>>"
    graph.add_node(pydot.Node(str(self.s.__hash__()),label=label,shape='plaintext'))

class MXOperationArtist(DotArtist):
  def draw(self):
    k = self.s
    graph = self.graph
    dep = [k.getDep(i) for i in range(k.getNdeps())]
    if not(k.isCommutative()):
      # Non-commutative operators are represented by 'record' shapes.
      # The dependencies have different 'ports' where arrows should arrive.
      s = getOperatorRepresentation(self.s,["| <f0> | ", " | <f1> |"])
      if s.startswith("(|") and s.endswith("|)"):
        s=s[2:-2]
      
      graph.add_node(pydot.Node(str(k.__hash__()),label=s,shape='Mrecord'))
      for i,n in enumerate(dep):
        graph.add_edge(pydot.Edge(str(n.__hash__()),str(k.__hash__())+":f%d" % i))
    else: 
     # Commutative operators can be represented more compactly as 'oval' shapes.
      s = getOperatorRepresentation(k,[".", "."])
      if s.startswith("(.") and s.endswith(".)"):
        s=s[2:-2]
      if s.startswith("(") and s.endswith(")"):
        s=s[1:-1]
      self.graph.add_node(pydot.Node(str(k.__hash__()),label=s,shape='oval'))
      for i,n in enumerate(dep):
        self.graph.add_edge(pydot.Edge(str(n.__hash__()),str(k.__hash__())))

class MXMultiplicationArtist(DotArtist):
  def draw(self):
    k = self.s
    graph = self.graph
    dep = [k.getDep(i) for i in range(k.getNdeps())]

    # Non-commutative operators are represented by 'record' shapes.
    # The dependencies have different 'ports' where arrows should arrive.
    s = "mul(| <f0> | , | <f1> | )"

    graph.add_node(pydot.Node(str(k.__hash__()),label=s,shape='Mrecord'))
    for i,n in enumerate(dep):
      graph.add_edge(pydot.Edge(str(n.__hash__()),str(k.__hash__())+":f%d" % i))
        
class SXMatrixArtist(DotArtist):
  def draw(self):
    s = self.s
    graph = self.graph
    sp = s.sparsity()
    row = sp.getRow()
      
    # The Matrix grid is represented by a html table with 'ports'
    label = '<<TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0">'
    for i in range(s.size1()):
      label+="<TR>"
      for j in range(s.size2()):
        k = sp.getNZ_const(i,j)
        if k==-1:
          label+="<TD>.</TD>"
        else:
          sx = s.at(k)
          if self.shouldEmbed(sx):
            label+="<TD BGCOLOR='#eeeeee'>%s</TD>" % str(sx)
          else:
            self.graph.add_edge(pydot.Edge(str(sx.__hash__()),"%s:f%d%d" % (str(self.s.__hash__()), i,j)))
            label+="<TD PORT='f%d%d' BGCOLOR='#eeeeee'> <font color='#666666'>(%d,%d)</font> </TD>" % (i,j,i,j)
      label+="</TR>"
    label+="</TABLE>>"
    graph.add_node(pydot.Node(str(self.s.__hash__()),label=label,shape='plaintext'))
    
  def shouldEmbed(self,sx):
    return len(self.invdep[sx]) == 1 and sx.isLeaf()
    
class SXLeafArtist(DotArtist):
  def draw(self):
    if len(self.invdep[self.s]) == 1:
      master = list(self.invdep[self.s])[0]
      if hasattr(self.artists[master],'shouldEmbed'):
        if self.artists[master].shouldEmbed(self.s):
          return
    style = "solid" # Symbolic nodes are represented box'es
    if self.s.isConstant():
      style = "bold" # Constants are represented by bold box'es
    self.graph.add_node(pydot.Node(str(self.s.__hash__()),label=str(self.s),shape="box",style=style)) 

class SXNonLeafArtist(DotArtist):
  def draw(self):
    k = self.s
    graph = self.graph
    dep = [k.getDep(i) for i in range(k.getNdeps())]
    if not(k.isCommutative()):
      # Non-commutative operators are represented by 'record' shapes.
      # The dependencies have different 'ports' where arrows should arrive.
      s = getOperatorRepresentation(self.s,["| <f0> | ", " | <f1> |"])
      if s.startswith("(|") and s.endswith("|)"):
        s=s[2:-2]
      
      graph.add_node(pydot.Node(str(k.__hash__()),label=s,shape='Mrecord'))
      for i,n in enumerate(dep):
        graph.add_edge(pydot.Edge(str(n.__hash__()),str(k.__hash__())+":f%d" % i))
    else: 
     # Commutative operators can be represented more compactly as 'oval' shapes.
      s = getOperatorRepresentation(k,[".", "."])
      if s.startswith("(.") and s.endswith(".)"):
        s=s[2:-2]
      if s.startswith("(") and s.endswith(")"):
        s=s[1:-1]
      self.graph.add_node(pydot.Node(str(k.__hash__()),label=s,shape='oval'))
      for i,n in enumerate(dep):
        self.graph.add_edge(pydot.Edge(str(n.__hash__()),str(k.__hash__())))
        
        
  
def createArtist(node,dep={},invdep={},graph=None,artists={}):
  if isinstance(node,SXMatrix):
    return SXMatrixArtist(node,dep=dep,invdep=invdep,graph=graph,artists=artists)
  elif isinstance(node,SX):
    if node.isLeaf():
      return SXLeafArtist(node,dep=dep,invdep=invdep,graph=graph,artists=artists)
    else:
      return SXNonLeafArtist(node,dep=dep,invdep=invdep,graph=graph,artists=artists)
  elif isinstance(node,MX):
    if node.isSymbolic():
      return MXSymbolicArtist(node,dep=dep,invdep=invdep,graph=graph,artists=artists)
    elif node.isBinary() or node.isUnary():
      return MXOperationArtist(node,dep=dep,invdep=invdep,graph=graph,artists=artists)
    elif node.isMultiplication():
      return MXMultiplicationArtist(node,dep=dep,invdep=invdep,graph=graph,artists=artists)
    elif node.isConstant():
      return MXConstantArtist(node,dep=dep,invdep=invdep,graph=graph,artists=artists)
  else:
    raise Exception("Cannot create artist for %s" % str(type(s)))
        
def dotgraph(s,direction="BT"):
  """
  Creates and returns a pydot graph structure that represents an SX or SXMatrix.
  
  direction   one of "BT", "LR", "TB", "RL"
  """
  
  # Get the dependencies and inverse dependencies in a dict
  dep, invdep = dependencyGraph(s,{},{})
  
  allnodes = set(dep.keys()).union(set(invdep.keys()))
  
  print "a", set(dep.keys()), [i.__hash__() for i in dep.keys()]
  print "b", set(invdep.keys()), [i.__hash__() for i in invdep.keys()]
  #print "allnodes", allnodes, [i.__hash__() for i in allnodes]
  
  #return None
  
  artists = {}
  
  graph = pydot.Dot('G', graph_type='digraph',rankdir=direction)
    
  for node in allnodes:
    artists[node] = createArtist(node,dep=dep,invdep=invdep,graph=graph,artists=artists)
    
  for artist in artists.itervalues():
    if artist is None: continue
    artist.draw()
  
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
