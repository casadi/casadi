from casadi import *
import copy

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
      
# fun must a function mapping from (flatcounter, item) to item    
def map_nested_list(fun,iterable):
  if isinstance(iterable,list): 
    ret = copy.deepcopy(iterable)
    for j in iter_flatten_mutator(ret,fun):
      pass
    return ret
  else:
    return fun(0,iterable)
 
class Variables(object):
    """
    A simple to use Variables container class
    """
    def __init__(self):
        self._d = dict()
        self._d_ = dict()
        self._orderflag = True
        self._offset = 0
        self._type = "SX"

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
        if name.startswith('o_') and len(name)>2:
            return self.getOffset(name[2:])
        if name.startswith('I_') and len(name)>2:
            return self.getIndex(name[2:])
        if name.startswith('i_') and len(name)>2:
            return self.getindex(name[2:])
        if name.startswith('iv_') and len(name)>2:
            return list(self.getindex(name[2:]))
        if name.endswith('_'):
            if isinstance(self._d[name[:-1]],Variables):
                raise Exception("Variable %s has no numerical value, because it is a Variables object itself." % name[:-1])
            self.doUpdates_(name[:-1])
            return self._d_[name[:-1]]
        if name in self._d:
          return self._d[name]
        else:
          raise Exception("Variable " + name + " does not exist. Existing ones are: " + ",".join(self.getOrder()))
           
    def __setattr__(self,name,value):
        """
        have a leading underscore to bypass the dictionary
        
        f.foo_ = 
        
        """
        if name.startswith('_'):
            object.__setattr__(self,name,value)
            return
        if name.endswith('_'): 
            getattr(self,name).set(value)
            return
        self._d[name] = value
        if isinstance(value,MX) or self._type == "MX":
           self._type = "MX"
           self.createParent()
        
        if not(isinstance(value,Variables)):
            self._d_[name] = self.get_(value)
            
    def createParent(self):
        self._V = msym("V[" + ",".join(self.getOrder()) + "]",self.getSize())
        for k in self.getOrder():
            obj = self._d[k]
            if isinstance(obj,Variables):
                raise Exception("Not implemented")
            i = flatten(getattr(self,'i_'+k))
            self._d[k] = map_nested_list(lambda c,x: self._V[i[c]],obj)
   
            
    def get_(self,value):
      return map_nested_list(lambda c,x: DMatrix(self.getSparsity(x),0),value)
            
    def doUpdates_(self,name):
      if isinstance(self._d[name],list):
        if not(len(self._d[name])==len(self._d_[name])):
          self._d_[name] = self.get_(self._d[name])
      
    def getindex(self,name):
        if not(name in self._d):
            raise Exception(("Variable %s not found. " % name)  + "\n" + "Available variables are: " + str(self.getOrder()))
            
        offsets = flatten(self.getOffset(name))
        return  map_nested_list(lambda c,x: self.getIMatrix(x,offsets[c]),self._d[name])
            
    def getIMatrix(self,obj,offset=0):
      s = self.getSize(obj)
      
      i = IMatrix(self.getSparsity(obj),range(s))
      ivec = vecNZ(i)
      i[ivec] = IMatrix(range(offset,s+offset))
      return i
            
    def getIndex(self,name):
        if not(name in self._d):
            raise Exception(("Variable %s not found. " % name)  + "\n" + "Available variables are: " + str(self.getOrder()))
            
        i = 0
        for k in self.getOrder():
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
            raise Exception(("Variable %s not found. " % name)  + "\n" + "Available variables are: " + str(self.getOrder()))
        i = self._offset

        for k in self.getOrder():
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
            return 0
                
    def getOrder(self):
        """
        Returns a list of variable names.
        This list defines the order of variables.
        """
        return sorted(self._d.iterkeys())
        
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
        if self._type == "MX":
          return self._V
          
        l = []
        for k in self.getOrder():
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
        if self._type == "MX":
          return self._V
        l = []
        for k in self.getOrder():
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
        l = []
        for k in self.getOrder():
            obj = self._d[k]
            if isinstance(obj,Variables):
                l.append(obj.veccat_())
            elif isinstance(obj,list):
                self.doUpdates_(k)
                l.append(vecNZcat(self._d_[k]))
            else:
                l.append(self._d_[k])
        return vecNZcat(flatten(l))
        
    def veccat_(self):
        """
        Returns the result of a veccat operation to the list of all variables values
        """
        l = []
        for k in self.getOrder():
            obj = self._d[k]
            if isinstance(obj,Variables):
                l.append(obj.veccat_())
            elif isinstance(obj,list):
                self.doUpdates_(k)
                l+=self._d_[k]
            else:
                l.append(self._d_[k])
        return veccat(flatten(l))
        
    def __str__(self):
        keys = self.getOrder()
        s=''
        s+= "Container holding %d variables.\n" % len(keys)
        for i,k in enumerate(keys):
            s+= ("%2d. " % i ) + k + " = " +  str(self._d[k]) + "\n"
        return s
   
    @property
    def shape(self):
       return (self.getNumel(),1)
       
       
