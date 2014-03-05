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


from casadi import *
import types

def frozen(b):
  def frozen_dec(method):
    def wrapper(self,*args,**kwargs):
      if self._frozen!=b:
        raise Exception("This operation can only be performed on %s instance." % ("a frozen" if b else "an unfrozen"))
      else:
        return method(self,*args,**kwargs)
    return wrapper
  return frozen_dec

class xcollection():
  @property
  def shape(self):
    return self[...].shape
    
  def cat(self):
    if not hasattr(self,"catted"):
      catter = veccat if self._keepZeros else vecNZcat
      self.catted = self._postcatmodifier(catter([self._modifier(n) for n in self._tree.traverse()]))
    return self.catted 
      
  def enumerate(self):
    for n in self._tree.traverse():
      yield (n.h,self._modifier(n))
      
  def items(self):
    for _,v in self.enumerate():
      yield v

class xdict(dict,xcollection):
  def __getitem__(self,name):
    if name is Ellipsis:
      return self.cat()
    elif isinstance(name,tuple):
      if len(name)==1:
        return self.__getitem__(name[0])
    else:
      if not name in self:
        raise Exception("Key '%s' not found. Only these are available: %s" % (name,str(self.keys())))
      return dict.__getitem__(self,name)
 
  def __getattr__(self,key):
    if key in self.__dict__:
      return object.__getattribute__(self,name)
    else:
      return self[key]

class xlist(list,xcollection): 
  def __getitem__(self,name):
    if name is Ellipsis:
      return self.cat()
    else:
      return list.__getitem__(self,name)
      
      

class Node:
  def __init__(self):
    self.h = '.'
    self.children = []
    self.super = None
    self._frozen = False
    self.i = 0
    self.exhausted = False
    self.labels = None
  
  @frozen(True)
  def copy(self):
    return self
    
  @frozen(False)
  def add(self,c):
    if isinstance(c,list):
      self.children+=c
      for i in c:
        i.parent = self
    else:
      self.children+=[c]
      c.parent = self
      
  def __str__(self):
    if len(self.children)==0:
      return str(self.h)+("X" if self.exhausted else "")
    else:
      return "(" + ",".join([str(c) for c in self.children]) + "|" + str(self.i) + ":" + ("X" if self.exhausted else "") + ("S" if self.super is not None else "")  +  ")"
  
  def freeze(self,h=[]):
    self._frozen = True
    self.h = h
    for i,v in enumerate(self.children):
      h_ = h
      labels = self.labels if self.labels is not None else range(len(self.children))
      if labels[i] is not None:
        h_ = h_ + [labels[i]]
      v.freeze(h_)
    self.reset()
      
  def traverseLeafs(self,maxdepth=1000):
    if maxdepth<=0:
      return
    if len(self.children)==0:
      yield self
    else:
      for v in self.children:
        for c in v.traverseLeafs(maxdepth-1):
          yield c
          
  def getGeneration(self,parent):
    gc = 0
    this = self
    while this is not parent:
      if not hasattr(this,'parent'):
        raise Exception("getGeneration: supplied node is not a parent")
      this = this.parent
      gc +=1
    return gc
          
  def traverseNodes(self,maxdepth=1000):
    if maxdepth<0:
      return
    if len(self.children)>0:
      yield self
      for v in self.children:
        for c in v.traverseNodes(maxdepth-1):
          yield c
             
  def reset(self):
    """
      Reset all counters
    """
    for v in self.traverseLeafs():
      v.exhausted = False
      
    for v in self.traverseNodes():
      v.i = 0
      v.exhausted = False
      
  @frozen(True)
  def structure(self,modifier = lambda x:x,decorate = lambda t,x:x,postcatmodifier=lambda x:x):
    if len(self.children)==0:
      return modifier(self)
    d = None
    if self.labels is not None:
      d = xdict()
      d._modifier = modifier
      d._postcatmodifier = postcatmodifier

      for i,v in enumerate(self.children):
        if self.labels[i] is None:
          d.update(v.structure(modifier=modifier,decorate=decorate,postcatmodifier=postcatmodifier))
        else:
          d[self.labels[i]] = v.structure(modifier=modifier,decorate=decorate,postcatmodifier=postcatmodifier)
      d = decorate(self,d)
    else:
      d = xlist()
      d._modifier = modifier
      d._postcatmodifier = postcatmodifier
      for v in self.children:
          d.append(v.structure(modifier=modifier,decorate = decorate,postcatmodifier=postcatmodifier))
      d = decorate(self,d)
    return d
    
  @frozen(True)
  def traverse(self):
    """
      Does smart tree traversal
    """
    # Reset all counters
    self.reset()
    
    # Locate the very first leaf (depth first)
    leaf = self
    while len(leaf.children)>0:
      leaf = leaf.children[leaf.i]

    # As long as we have a non-None leaf, yield it and look up the next leaf
    while leaf:
      # Yield the leaf and mark it as exhausted
      yield leaf
      leaf.exhausted = True
      
      # First stage - bubble up
      this = leaf.parent
      while 1:
        # Update counter and re-calculate exhaustedness
        this.i = (this.i + 1) % len(this.children)
        this.exhausted = all([i.exhausted for i in this.children]) 
        if this.super is not None: # Follow super routes
          this = this.super
        elif this.exhausted:       # Bubble up if exhausted
          if hasattr(this,'parent'):
            this = this.parent
          else: 
            return # No more bubbling up possible; done yielding
        else:
          break # The current node will serve as upper node

      # Second stage, go down   
      while len(this.children)>0:
        # Select the first child
        this = this.children[this.i]
        while this.exhausted: # Bubble up if exhausted
          if not hasattr(this,'parent'):
            return # No more bubbling up possible; done yielding
          this = this.parent
          this.i = (this.i + 1) % len(this.children)
          this.exhausted = all([i.exhausted for i in this.children])
        
      # Remember found leaf for next iteration
      leaf = this
      
class Mapping:
  def __init__(self):
    self.data = {}
    self.nameOrder = []
    self.orderInstructions = {}
    
  def add(self,name,data,metadata = None):
    if name in self.data:
      raise Exception("Name '%s' is already used." % name)
    else:
      self.data[name] = data
      self.nameOrder += [name]
  
  def freeze(self):
    if not hasattr(self,'order'):
      self.order = self.nameOrder
    return FrozenMapping(self)
    
  def setOrder(self,order):
    ordernames = []
    for v in order:
      if isinstance(v,str):
        ordernames.append(v)
      elif isinstance(v,tuple):
        for e in v:
          if isinstance(e,str):
            ordernames.append(e)
          elif isinstance(e,tuple):
            ordernames.append(e[0])
      else:
        raise Exception("setOrder error: %s" % v)
    for i in ordernames:
      if not i in self.nameOrder:
        raise Exception("setOrder(order): '%s' is not found in data dict. Available fields are: %s" % (str(i),str(self.nameOrder)))
    for i in self.nameOrder:
      if not i in ordernames:
        raise Exception("setOrder(order): '%s' from data dict was not found in supplied fields. Supplied fields are: %s" % (str(i),str(ordernames)))
    
    self.order = order

def getLeafSize(node):
  return node.size()
    
class FrozenMapping:
  def __init__(self,mapping):
    self.data = mapping.data
    self.order = mapping.order
    self.orderInstructions = mapping.orderInstructions
    self.tree = Node()
    self.tree.labels = []
    for k in self.order:
      if isinstance(k,str):
        self.createNodes(mapping.data[k],self.tree)
        self.tree.labels.append(k)
      elif isinstance(k,tuple):
        if all([isinstance(e,str) or isinstance(e,tuple) for e in k]):
          n = Node()
          n.labels = []
          self.tree.add(n)
          self.tree.labels.append(None)
          for e in k:
            name = e if isinstance(e,str) else e[0]
            axis = 0
            if isinstance(e,tuple):
              for j in e:
                if isinstance(j,int):
                  axis = j
            n.labels.append(name)
            self.createNodes(mapping.data[name],n)
          for node in n.traverseNodes(maxdepth=axis+1):
            if node.getGeneration(n)==axis+1:              
              node.super = n
        else:
          name = k[0]
          self.createNodes(mapping.data[name],self.tree)
          self.tree.labels.append(name)
      else:
        raise Exception("panic")

    self.tree.freeze()
    
    offset = 0
    cnt = 0
    for n in self.tree.traverse():
      n.index = IMatrix(n.data.sparsity(),range(offset,offset+getLeafSize(n.data)))
      n.offset = offset
      n.count = cnt
      cnt += 1
      offset += getLeafSize(n.data)
    self.size = offset
    self.count = cnt
    self.datacount = len(self.data)
      
  def createNodes(self,data,parent):
    if isinstance(data,list):
      n = Node()
      parent.add(n)
      for e in data:
        self.createNodes(e,n)
    elif isinstance(data,Collection):
      parent.add(data._mapping.tree.copy()) 
    elif isinstance(data,FrozenMapping):
      parent.add(data.tree.copy())
    else:
      if isinstance(data,SXElement):
        data = SX(data)
      n = Node()
      n.data = data
      parent.add(n)
      
class Mapper:
  def __init__(self,mapping,modifier=None,decorate=None,postcatmodifier=None):
    self.mapping = mapping
    self.structure = mapping.tree.structure(modifier=modifier,decorate=decorate,postcatmodifier=postcatmodifier)

  def __getitem__(self, name):
    return self.structure.__getitem__(name)
  
class Collection(object):
    def __init__(self,keepZeros=True):
      self._frozen = False
      self._keepZeros = keepZeros
      self._mapping = Mapping()

    @frozen(False)
    def freeze(self):
      self._frozen = True
      self._mapping = self._mapping.freeze()
      

      def catter(tree,iter):
        iter._tree=tree
        iter._keepZeros = self._keepZeros
        return iter
        
      self._i  = Mapper(self._mapping,lambda x: x.index,decorate=catter,postcatmodifier=lambda x:x)
      self._iv = Mapper(self._mapping,lambda x: list(x.index),decorate=catter,postcatmodifier=lambda x:list(x))
      self._o = Mapper(self._mapping,lambda x: x.offset,decorate=catter,postcatmodifier=lambda x:x)
      self._d = Mapper(self._mapping,lambda x: x.data,decorate=catter,postcatmodifier=lambda x:x)
      
    def __getattr__(self,name): 
      """
        f.foo   | f["foo"]    returns the 'foo' entry in the variable dictionary
        
        f.i_foo | f.i_["foo"] returns the IMatrix index 
        f.iv_foo same as list(f.i_foo)
        f.I_foo returns the flat index
        f.o_foo returns the offset 
        f._foo  returns the foo attribute of the object f
        f.foo_  returns a DMatrix with the same sparsity as the 'foo entry'
        
      """
      if name in self.__dict__:
          return object.__getattribute__(self,name)
      if not(self._frozen):
        raise Exception("You must first freeze() this Collection instance before you can access its members ('%s' in this case). This is for your own protection" % name)
      if name.startswith('o_') and len(name)>2:
          return object.__getattribute__(self,'_'+name[:1])[name[2:]]
      if name.startswith('I_') and len(name)>2:
          return object.__getattribute__(self,'_'+name[:1])[name[2:]]
      if name.startswith('i_') and len(name)>2:
          return object.__getattribute__(self,'_'+name[:1])[name[2:]]
      if name.startswith('iv_') and len(name)>3:
          return object.__getattribute__(self,'_'+name[:2])[name[3:]]

      return self._d[name]

    def __setattr__(self,name,value):
        """
        have a leading underscore to bypass the dictionary
        
        f.foo_ = 
        
        """
        if name.startswith('_'):
            object.__setattr__(self,name,value)
            return
        if self._frozen:
          raise Exception("This Collection instance is frozen. You cannot add more members ('%s' in this case). This is for your own protection." % name)
        self[name] = value 
           
    @frozen(True)
    def __getitem__(self, name):
      return self._d[name]

    @frozen(False)
    def __setitem__(self, name, value):
      self._mapping.add(name,value)
         
    @property
    @frozen(True)
    def shape(self):
       return (self._mapping.size,1)

    @property
    @frozen(True)
    def size(self):
       return self._mapping.size

    def __str__(self):
        if self._frozen:
          s=''
          s+= "Frozen container holding %d variables.\n" % self._mapping.datacount
          s+="Order: %s\n" % str(self._mapping.order)
        else:
          s=''
          s+= "Unfrozen container holding %d variables.\n" % self._mapping.datacount
        for k in self._mapping.data.iterkeys():
          s+= k + " = " +  str(self._mapping.data[k]) + "\n"
        return s
       
    @frozen(False)
    def setOrder(self,order):
      self._mapping.setOrder(order)
      
    @frozen(True)
    def numbers(self,init=0):
      return DMatrix(self.size,1,init)
      
    @frozen(False)
    def cat(self):
      return self._d[...]

