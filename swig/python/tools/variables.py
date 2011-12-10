from casadi import *

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
        if name.endswith('_'):
            if isinstance(self._d[name[:-1]],Variables):
                raise Exception("Variable %s has no numerical value, because it is a Variables object itself." % name[:-1])
            self.doUpdates_(name[:-1])
            return self._d_[name[:-1]]
        return self._d[name]
           
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
        self._V = msym("V",self.getSize())
        for k in self.getOrder():
            obj = self._d[k]
            if isinstance(obj,Variables):
                raise Exception("Not implemented")
            elif isinstance(obj,list):
                for i in range(len(obj)):
                    self._d[k] = V[getattr(self,'i_'+k)[i]]
            else:
                print getattr(self,'i_'+k)
                print self._V.shape, self._V.size()
                print self._V[getattr(self,'i_'+k)] 
                
                self._d[k] = self._V[getattr(self,'i_'+k)]  
                print "foo"      
            
    def get_(self,value):
        if isinstance(value,list):
            ret = []
            for v in value:
              ret.append(self.get_(v))
            return ret
        else:
            return DMatrix(self.getSparsity(value),0)
            
    def doUpdates_(self,name):
      if isinstance(self._d[name],list):
        if not(len(self._d[name])==len(self._d_[name])):
          self._d_[name] = self.get_(self._d[name])
      
    def getindex(self,name):
        if not(name in self._d):
            raise Exception(("Variable %s not found. " % name)  + "\n" + "Available variables are: " + str(self.getOrder()))
            
        if isinstance(self._d[name],list):
            offset = self.getOffset(name)[0]
            ret = []
            for v in self._d[name]:
              s = self.getSize(v)
              i = IMatrix(self.getSparsity(v),range(s))
              ivec = vecNZ(i)
              i[ivec] = IMatrix(range(offset,s+offset))
              ret.append(i)
              offset += s
            obj = self._d[name]
            
            return ret
        else:
            offset = self.getOffset(name)
            obj = self._d[name]
            s = self.getSize(obj)
            
            i = IMatrix(self.getSparsity(obj),range(s))
            ivec = vecNZ(i)
            i[ivec] = IMatrix(range(offset,s+offset))
            return i
            
    def getIndex(self,name):
        if not(name in self._d):
            raise Exception(("Variable %s not found. " % name)  + "\n" + "Available variables are: " + str(self.getOrder()))
            
        if isinstance(self._d[name],list):
            offset = self.getOrder().index(name)
            return range(offset,offset+len(self._d[name]))
        else:
            i = 0
            for k in self.getOrder():
               if k is name:
                  break
               if isinstance(self._d[k],list):
                  i+=len(self._d[k])
               else:
                  i+=1

            return i
             
    def setOffset(self,value):
        self._offset = value
        
    def getOffset(self,name):
        if not(name in self._d):
            raise Exception(("Variable %s not found. " % name)  + "\n" + "Available variables are: " + str(self.getOrder()))
        i = self._offset
        if isinstance(self._d[name],list):
          for k in self.getOrder():
              if k == name:
                  ret = []
                  for v in self._d[name]:
                    ret.append(i)
                    i += self.getSize(v)
                  return ret
              else:
                  i += self.getSize(self._d[k])
        else:
          for k in self.getOrder():
              if k == name:
                  return i
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
        
    def getNumel(self):
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
        return veccat(l)
        
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
        return vecNZcat(l)

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
        return vecNZcat(l)
        
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
        return veccat(l)
        
    def __str__(self):
        keys = self.getOrder()
        s=''
        s+= "Container holding %d variables.\n" % len(keys)
        for i,k in enumerate(keys):
            s+= ("%2d. " % i ) + k + " = " +  str(self._d[k]) + "\n"
        return s
       
