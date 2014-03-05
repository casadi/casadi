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
import warnings
from casadi import *
import copy
from types import *

def flatten(l):
  return [i for i in iter_flatten(l)]
 
def iter_flatten(iterable):
  if isinstance(iterable,list):
    for e in iterable:
      if isinstance(e, (list, tuple)):
        for f in iter_flatten(e):
          yield f
      else:
        yield e
  else:
    yield iterable
    
def iter_flatten_list(iterable):
  if isinstance(iterable,list):
    for e in iterable:
      if isinstance(e, list):
        for f in iter_flatten(e):
          yield f
      else:
        yield e
  else:
    yield iterable
    
def iter_flatten_mutator(iterable,mutator=lambda i,x:x,counter=0):
  if isinstance(iterable,list):    
    for i in range(len(iterable)):
      e = iterable[i] 
      if isinstance(e, (list, tuple)):
        for c,f in iter_flatten_mutator(e,mutator,counter=counter):
          counter = c
          yield (counter,f)
          counter+=1
      else:
        iterable[i] = mutator(counter,e)
        yield (counter,e)
        counter+=1
  else:
     yield (counter,iterable)
     counter+=1 
      
def iter_flatten_mutator_hierarchical(iterable,mutator=lambda i,j,x:x,counter=0,hierarchy=()):
  if isinstance(iterable,list):    
    for i in range(len(iterable)):
      h = hierarchy + (i,)
      e = iterable[i] 
      if isinstance(e, (list, tuple)):
        for c,f in iter_flatten_mutator_hierarchical(e,mutator,counter=counter,hierarchy=(i,) + hierarchy):
          counter = c
          yield (counter,f)
          counter+=1
      else:
        iterable[i] = mutator(counter,h,e)
        yield (counter,e)
        counter+=1
  else:
     yield (counter,iterable)
     counter+=1 
      
# fun must a function mapping from (flatcounter, item) to item    
def map_nested_list(fun,iterable):
  if isinstance(iterable,list): 
    ret = copy.deepcopy(iterable)
    for j in iter_flatten_mutator(ret,fun):
      pass
    return ret
  else:
    return fun(0,iterable)
    
# fun must a function mapping from (flatcounter, hierarchicalcountertuple, item) to item    
def map_nested_list_hierarchical(fun,iterable):
  if isinstance(iterable,list): 
    ret = copy.deepcopy(iterable)
    for j in iter_flatten_mutator_hierarchical(ret,fun):
      pass
    return ret
  else:
    return fun(0,iterable)
 
class Variables(object):
    """
    A simple to use Variables container class
    """
    def __init__(self):
        warnings.warn("Variables class is deprecated. Use Collection instead", DeprecationWarning)
        self._d = dict()
        self._d_ = dict()
        self._orderflag = True
        self._offset = 0
        self._type = "SX"
        self._frozen = False
        
    def unfrozencopy(self):
      c = Variables()
      c._offset = copy.copy(self._offset)
      c._d = copy.copy(self._d)
      c._d_ = copy.copy(self._d_)
      c._type = copy.copy(self._type)
      c._orderflag = copy.copy(self._orderflag)
      return c
        
    def freeze(self,parent=True,bare=False):
      self._createParent = parent
      self._frozen = True
      self._order = sorted(self._d.iterkeys())
      for k in self._order:
        if isinstance(self._d[k],Variables):
          self._d[k].freeze()
      if self._type == "MX" and self._createParent:
        self.createParent()
      if not(bare):
        self._numbers = Numbers(self,recycle=True)
        self.buildlookuptable()
                
    def buildlookuptable(self):
      self._reverselookup  = [None]*self.getSize()
      self._reverselookup2 = [None]*self.getSize()
      for k in self._order:
          obj = self._d[k]
          offset = self.getOffset(k)
          size = self.getSize(obj)
          if isinstance(obj,Variables):
            self._reverselookup[offset:offset+size]   = [ (k,)+i for i in obj._reverselookup]
            self._reverselookup2[offset:offset+size] = [ (k,)+i for i in obj._reverselookup2]
          elif isinstance(obj,list):
            result = map_nested_list_hierarchical(lambda flati, hieri, item: (hieri, item ) ,self.getindex(k))
            for hierarchy, imatrix in iter_flatten_list(result):
              for i,j in enumerate(list(imatrix)):
                self._reverselookup[j] = (k,)+ hierarchy + ((i,),)
                self._reverselookup2[j] = (k,)+ hierarchy + ((i,),)
          else:
           for i,j in enumerate(list(self.getindex(k))):
             self._reverselookup[j]  = (k,(i,))
             self._reverselookup2[j] = (k,(i,))
                
                
    def reverselookup(self,index):
      return self._reverselookup[index]
      
    @property
    def reverselookuptable(self):
      return self._reverselookup
      
    def lookup(self,index,obj = None):
      """
      
      index 
      
      ('x',)
      ('x',(1,))
      ('y',2,(1,))
      
      """
      if obj is None:
        obj = self
      
      if not(isinstance(index,tuple)):
        index = (index,)
      
      if len(index)==0:
        return obj
        
      first = index[0]
      rest = index[1:]
      
      if isinstance(obj,list):
        return self.lookup(rest,obj=obj[first]) 
      elif isinstance(obj,Variables):
        return self.lookup(rest,obj=obj.__getattr__(first)) 
      elif isinstance(obj,SXElement):
        if len(first)==1 and first[0]==0:
          return obj
        else:
          raise Exception("Invalid index")
      else:
        return self.lookup(rest,obj=obj.__getitem__(first)) 
      
    
    def hierarchicalIndexRepr(self,index,matrixstyle=False):
      s = ""
      for i in index:
        if isinstance(i,StringType):
          s+=i
        elif isinstance(i,TupleType):
          s+=str(list(i))
        elif isinstance(i,int):
          s+="[%d]" % i
      return s
      
    def getLabels(self):
      return [self.hierarchicalIndexRepr(i) for i in self._reverselookup2]
      
    def __getattr__(self,name):
        """
          f.foo   returns the 'foo' entry in the variable dictionary 
          f.i_foo returns the IMatrix index 
          f.iv_foo same as list(f.i_foo)
          f.I_foo returns the flat index
          f.o_foo returns the offset 
          f._foo  returns the foo attribute of the object f
          f.foo_  returns a DMatrix with the same sparsity as the 'foo entry'
          
        """
        if name in self.__dict__:
            return object.__getattribute__(self,name)
        if not(self._frozen):
          raise Exception("You must first freeze() this Variables instance before you can access its members ('%s' in this case). This is for your own protection" % name)
        if name.startswith('o_') and len(name)>2:
            return self.getOffset(name[2:])
        if name.startswith('I_') and len(name)>2:
            return self.getIndex(name[2:])
        if name.startswith('i_') and len(name)>2:
            return self.getindex(name[2:])
        if name.startswith('iv_') and len(name)>2:
            return list(self.getindex(name[3:]))
        if name.endswith('_'):
            if isinstance(self._d[name[:-1]],Variables):
                raise Exception("Variable %s has no numerical value, because it is a Variables object itself." % name[:-1])
            return getattr(self._numbers,name[:-1])
        if name in self._d:
          return self._d[name]
        else:
          raise Exception("Variable " + name + " does not exist. Existing ones are: " + ",".join(self._order))
           
    def __setattr__(self,name,value):
        """
        have a leading underscore to bypass the dictionary
        
        f.foo_ = 
        
        """
        if name.startswith('_'):
            object.__setattr__(self,name,value)
            return
        if name.endswith('_'): 
            getattr(self._numbers,name[:-1]).set(value)
            return
        if self._frozen:
          raise Exception("This Variables instance is frozen. You cannot add more members ('%s' in this case). This is for your own protection." % name)
        if isinstance(value,Variables):
          value = value.unfrozencopy()
          value.freeze()
        self._d[name] = value
        if (isinstance(value,MX) and value.isSymbolic()) or self._type == "MX":
           self._type = "MX"
        
            
    def createParent(self):
        self._V = MX.sym("V[" + ",".join(self._order) + "]",self.getSize())
        for k in self._order:
            obj = self._d[k]
            if isinstance(obj,Variables):
                raise Exception("Not implemented")
            i = flatten(getattr(self,'i_'+k))
            self._d[k] = map_nested_list(lambda c,x: self._V[i[c]],obj)
      
    def getindex(self,name):
        if not(name in self._d):
            raise Exception(("Variable %s not found. " % name)  + "\n" + "Available variables are: " + str(self._order))
            
        if not(hasattr(self,'_indexCache')):
          self._indexCache={}
          
        if name in self._indexCache:
          return self._indexCache[name]

        offsets = flatten(self.getOffset(name))
          
        if isinstance(self._d[name],Variables):
          res = IMatrix(range(offsets[0],offsets[0]+self._d[name].shape[0]))
        else:
          res = map_nested_list(lambda c,x: self.getIMatrix(x,offsets[c]),self._d[name])
          
        
        self._indexCache[name] = res
        
        return res
            
    def getIMatrix(self,obj,offset=0):
        
      s = self.getSize(obj)
      
      i = IMatrix(self.getSparsity(obj),range(s))
      ivec = vecNZ(i)
      i[ivec] = IMatrix(range(offset,s+offset))
      
      
      return i
            
    def getIndex(self,name):
        if not(name in self._d):
            raise Exception(("Variable %s not found. " % name)  + "\n" + "Available variables are: " + str(self._order))
            
        i = 0
        for k in self._order:
           if k is name:
              break
           if isinstance(self._d[k],list):
              i+=len(flatten(self._d[k]))
           else:
              i+=1 

        return map_nested_list(lambda c,x: c+i,self._d[name])
             
    def setOffset(self,value):
        self._offset = value
        
    def getOffset(self,name):
        if not(name in self._d):
            raise Exception(("Variable %s not found. " % name)  + "\n" + "Available variables are: " + str(self._order))
        i = self._offset

        for k in self._order:
            if k == name:
                sizes = [self.getSize(v) for v in iter_flatten(self._d[name])]
                offsets = [i] * (len(sizes)+1)
                
                for j in range(len(sizes)):
                  offsets[j+1] = offsets[j]+sizes[j]
                
                return map_nested_list(lambda c,x: offsets[c],self._d[name])
            else:
                i += self.getSize(self._d[k])

        
    def getSize(self,ob = None):
        if isinstance(ob,Variables):
            return ob.getSize()
        elif isinstance(ob,list):
            i = 0
            for v in ob:
                i+=self.getSize(v)
            return i
        elif hasattr(ob,'size'):
            return ob.size()
        elif ob is None:
            i = 0
            for v in self._d.itervalues():
                i+=self.getSize(v)
            return i
        else:
            return 1
        
    def getNumel(self,ob = None):
        if isinstance(ob,Variables):
            return ob.getNumel()
        elif hasattr(ob,'numel'):
            return ob.numel()
        elif isinstance(ob,list):
            i = 0
            for v in ob:
                i+=self.getNumel(v)
            return i
        elif ob is None:
            i = 0
            for v in self._d.itervalues():
                i+=self.getNumel(v)
            return i
        else:
            return 1
        
    def getSparsity(self, obj = None):
        if obj is None:
            return self.veccat().sparsity()
        elif hasattr(obj,'sparsity'):
            return obj.sparsity()
        else:
            return sp_dense(1,1)
        
    def veccat(self):
        """
        Returns the result of a veccat operation to the list of all variables
        """
        if not(self._frozen):
          raise Exception("You must first freeze() this Variables instance. This is for your own protection")
        if self._type == "MX" and self._createParent:
          return self._V
          
        l = []
        for k in self._order:
            obj = self._d[k]
            if isinstance(obj,Variables):
                l.append(obj.veccat())
            elif isinstance(obj,list):
                for v in obj:
                  l.append(v)
            else:
                l.append(self._d[k])
        return veccat(flatten(l))
        
    def vecNZcat(self):
        """
        Returns the result of a veccat operation to the list of all variables
        """
        if not(self._frozen):
          raise Exception("You must first freeze() this Variables instance. This is for your own protection")
        if self._type == "MX" and self._createParent:
          return self._V
        l = []
        for k in self._order:
            obj = self._d[k]
            if isinstance(obj,Variables):
                l.append(obj.vecNZcat())
            elif isinstance(obj,list):
                for v in obj:
                  l.append(v)
            else:
                l.append(self._d[k])
        return vecNZcat(flatten(l))

    def vecNZcat_(self):
        """
        Returns the result of a veccat operation to the list of all variables values
        """
        if not(self._frozen):
          raise Exception("You must first freeze() this Variables instance. This is for your own protection")
        return self._numbers.vecNZcat()
        
    def veccat_(self):
        """
        Returns the result of a veccat operation to the list of all variables values
        """
        if not(self._frozen):
          raise Exception("You must first freeze() this Variables instance. This is for your own protection")
        return self._numbers.veccat()
        
    def __str__(self):
        if self._frozen:
          keys = self._order
          s=''
          s+= "Frozen container holding %d variables.\n" % len(keys)
          for i,k in enumerate(keys):
              s+= ("%2d. " % i ) + k + " = " +  str(self._d[k]) + "\n"
          return s
        else:
          s=''
          s+= "Unfrozen container holding %d variables.\n" % len(keys)
          for k in self._d.iterkeys():
              s+= k + " = " +  str(self._d[k]) + "\n"
          return s
   
    @property
    def shape(self):
       return (self.getNumel(),1)
          
class Numbers(object):
  def __init__(self,variables=None, recycle = False):
    if (isinstance(variables,Numbers)):
      self._variables = variables._variables
      self._d_ = copy.deepcopy(variables._d_)
      self._numbers = copy.deepcopy(variables._numbers)
      return
      
    if (variables is None or not(isinstance(variables,Variables)) or not(variables._frozen)):
      raise Exception("You must supply the Numbers constructor with a Variables instance in frozen configuration.")
      
      
    self._variables = variables
    self._d_ = dict()
    self._numbers = dict()
    
    for k in self._variables._order:
      obj = self._variables._d[k]
      if isinstance(obj,Variables):
        if recycle:
          self._numbers[k] = obj._numbers
        else:
          self._numbers[k] = Numbers(obj)
      else:
        self._d_[k] = map_nested_list(lambda c,x: DMatrix(self._variables.getSparsity(x),0),obj)
      
  def __getattr__(self,name):
        """
          f.foo   returns the 'foo' DMatrix entry in the variable dictionary 
          
        """
        if name in self.__dict__:
            return object.__getattribute__(self,name)
            
        if name in self._numbers:
          return self._numbers[name]
            
        if name in self._d_:
          return self._d_[name]
        else:
          raise Exception("Variable " + name + " does not exist. Existing ones are: " + ",".join(self._variables._order))
           
  def __setattr__(self,name,value):
        """
        have a leading underscore to bypass the dictionary
        
        f.foo 
        
        """
        
        if name.startswith('_'):
            object.__setattr__(self,name,value)
            return
        
        getattr(self,name).set(value)
            
  def vecNZcat(self):
      """
      Returns the result of a veccat operation to the list of all variables values
      """
      l = []
      for k in self._variables._order:
          obj = self._variables._d[k]
          if isinstance(obj,Variables):
              l.append(self._numbers[k].veccat())
          elif isinstance(obj,list):
              l.append(vecNZcat(self._d_[k]))
          else:
              l.append(self._d_[k])
      return vecNZcat(flatten(l))
      
  def veccat(self):
      """
      Returns the result of a veccat operation to the list of all variables values
      """
      l = []
      for k in self._variables._order:
          obj = self._variables._d[k]
          if isinstance(obj,Variables):
              l.append(self._numbers[k].veccat())
          elif isinstance(obj,list):
              l+=self._d_[k]
          else:
              l.append(self._d_[k])
      return veccat(flatten(l))
      
  def setAll(self,value):
    for k in self._variables._order:
      getattr(self,k).setAll(value)
