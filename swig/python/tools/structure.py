from casadi import *

import operator

# StructIndex :tuple/list of strings
# canonicalIndex : tuple/list of string or numbers
# powerIndex: tuple/list of string, numbers, lists, slices, dicts

# flatIndex

# Primitive helpers
def lpack(L): return [[x] for x in L]
  
def combine(*args):
  if len(args)==0:
    return [[]]
  else:
    return [a + b for a in args[0] for b in combine(*args[1:])]

def listindices(dims,nest=False):
  if len(dims)==0:
    if nest:
      return [[[]]]
    else:
      return [[]]
  else:
    tail = listindices(dims[1:])
    if nest:
      return [combine([[i]],tail) for i in range(dims[0])]
    else:
      return combine(lpack(range(dims[0])),tail)

def intersperseIt(*args):
  iterators = map(iter,args)
  active = [True]*len(args)
  i = 0
  while any(active):
    try:
      yield iterators[i].next()
    except:
      active[i] = False
    i = (i + 1) % len(args)

def intersperse(*args):
   return list(intersperseIt(*args))
   
def canonicalIndexAncestors(ind):
  if len(ind)==0: return []
  isstring = lambda x: isinstance(x,str)
  return [ind] + canonicalIndexAncestors(ind[:-(map(isstring,ind[::-1]).index(True)+1)])

def canonical(ind,s):
  if ind < 0:
    return ind + s
  else:
    return ind
    
def flatten(e):
  if any(isinstance(i,list) for i in e):
    return sum(map(flatten,e),[])
  else:
    return e
    
# Decoraters

def properGetitem(f):
  """
    This decorator modifies a __getitem__/__setitem__ method such that it will always receive a tuple
  """
  def proper(self,mbt,*args):
    if not isinstance(mbt,tuple):
      mbt = (mbt,)
    return f(self,mbt,*args)
  return proper
  
# Enhanced standard classes
class SafeDict(dict):
  def __getitem__(self,k):
    if k not in self:
      raise Exception("Unknown keyword '%s'. Available entries: %s" % (k,str(self.keys())))
    return dict.__getitem__(self,k)
  
# Placeholder classes and instances
class Repeater:
  def __init__(self,e): self.e = e

def repeated(e):
  """
    From the arguemnt, constructs something that acts like a 'list' with the argument repeated the 'correct' number of times 
    
    s = struct_ssym([entry("x",repeat=6)])
    s["x"] = repeated(12)
    
  """
  return Repeater(e)
  
class NestedDictLiteral:
  """
    NestedDictLiteral will cause all dictionaries to become explicit recursively
  """
  
nesteddict = NestedDictLiteral()

# Casadi-independant Structure framework

def payloadUnpack(payload,i):
  if isinstance(i,str):
    raise Exception("Got string %s where number expected."% i)
  if isinstance(payload,list):
    if i>=len(payload):
      raise Exception("Rhs out of range. Got list index %s but rhs list is only of length %s." % (i,len(payload)))
    return payload[i]
  elif isinstance(payload,Repeater):
    return payload.e
  else:
    return payload

class StructEntry:
  def __init__(self,name,struct=None,data=None,dims=[]):
    self.name = name
    self.dims = dims
    self.struct = struct

  def __str__(self,compact=False):
     s=''
     if len(self.dims)>=1:
       s+= "repeated(%s): " % str(self.dims)
     if self.isPrimitive():
       s+=self.primitiveString()
     else:
       s+=self.struct.__str__(compact=True)
     return s
  
  __repr__ = __str__
     
  def primitiveString(self):
    return "*"
  
  def isPrimitive(self):
    return self.struct is None
    
  def traverseCanonicalIndex(self,nest=False,limit=1000):
    children = [[]] if (self.struct is None or limit==0) else self.struct.traverseCanonicalIndex(limit=limit-1)
    li = listindices(self.dims,nest)
    n = [[self.name]]
    if nest:
      return [combine(n,i,children) for i in li]
    else:
      return combine(n,li,children)
      
  def getStructEntryByStructIndex(self,structIndex):
    return self.struct.getStructEntryByStructIndex(structIndex)
      
  def traverseByPowerIndex(self,powerIndex,dims=None,canonicalIndex=(),dispatcher=None,payload=None):
    try:
      if dims is None: dims = self.dims
      # At the end of powerIndex, pending : are added automatically if dims is not exhausted
      if len(powerIndex)==0:  
        if len(dims)>0:
          return self.traverseByPowerIndex(
                   [slice(None) for i in dims],
                   dims=dims,
                   canonicalIndex=canonicalIndex,
                   dispatcher=dispatcher,
                   payload=payload
                 )
        else:
          return dispatcher(payload,canonicalIndex,entry=self)
          
      if len(dims)==0:
        if self.isPrimitive(): # Pass on remainder of powerIndex to dispatcher
          return dispatcher(payload,canonicalIndex,extraIndex=tuple(powerIndex),entry=self)
        else:
          return self.struct.traverseByPowerIndex(
                   powerIndex,
                   canonicalIndex=canonicalIndex,
                   dispatcher=dispatcher,
                   payload=payload
                 )
      else:
        p = powerIndex[0]
        s = dims[0]
        if isinstance(p,slice): # Expand slice
          p = range(*p.indices(s))
        if isinstance(p,int):
          return self.traverseByPowerIndex(
                   powerIndex[1:],
                   dims=dims[1:],
                   canonicalIndex=canonicalIndex+(canonical(p,s),),
                   dispatcher=dispatcher,
                   payload=payload
                 )
        elif isinstance(p,list):
          return [
                    self.traverseByPowerIndex(
                      powerIndex[1:],
                      dims=dims[1:],
                      canonicalIndex=canonicalIndex+(canonical(i,s),),
                      dispatcher=dispatcher,
                      payload = payloadUnpack(payload,i)
                    )
                 for i in p]
        elif isinstance(p,dict):
          raise Exception("powerIndex entry {} cannot be used in list context.")
        elif isinstance(p,NestedDictLiteral):
          return [
                    self.traverseByPowerIndex(
                      [p],
                      dims=dims[1:],
                      canonicalIndex=canonicalIndex+(canonical(i,s),),
                      dispatcher=dispatcher,
                      payload = payloadUnpack(payload,i)
                    )
                 for i in range(s)]
        elif callable(p):
          r = p(self.traverseByPowerIndex(
                powerIndex[1:],
                dims=dims,
                canonicalIndex=canonicalIndex,
                dispatcher=dispatcher.callableInner(),
                payload=payload
              ))
          return dispatcher.callableOuter(payload,canonicalIndex,extraIndex=None,entry=None,inner=r)
        else:
          raise Exception("I don't know what to do with this: %s" % str(p))
    except Exception as e:
      raise Exception("Error occured in entry context with powerIndex %s, at canonicalIndex %s:\n%s" % (str(powerIndex),str(canonicalIndex),str(e)))
  
class Structure:
  def __init__(self,entries,order=None):
    self.entries = entries
    
    self.order = [e.name for e in self.entries] if order is None else order
    self.keyslist = sum([ list(i) if isinstance(i,tuple) else list([i]) for i in self.order],[])
    
    self.dict = SafeDict([(e.name,e) for e in self.entries])
    
    for e in self.order:
      if isinstance(e,str):
        if e not in self.dict:
          raise Exception("Order '%s' is invalid." % e)
      elif isinstance(e,tuple):
        for ee in e:
          if ee not in self.dict:
            raise Exception("Order '%s' is invalid." % ee)
       
  def keys(self):
    return self.keyslist

  def __str__(self,compact=False):
     s=''
     if compact:
       s+= "{" + ",".join(k + ": " +  v.__str__(compact=True) for k,v in self.dict.items()) + "}"
     else:
       s+= "Structure holding %d entries.\n" % len(self.dict)
       s+="  Order: %s\n" % str(self.order)
       for k,v in self.dict.items():
          s+= "  " + k + " = " +  v.__str__(compact=True) + "\n"
     return s
     
  __repr__ = __str__
        
  def traverseCanonicalIndex(self,limit=1000):
    ret = []
    for d in self.order:
      if isinstance(d,tuple):
        for v in intersperse(*[self.dict[de].traverseCanonicalIndex(True,limit=limit-1) for de in d]):
          ret += v 
      else:
        ret += self.dict[d].traverseCanonicalIndex(limit=limit-1)
    return ret

  def getStructEntryByStructIndex(self,structIndex):
    e = self.dict[structIndex[0]]
    if len(structIndex)>1:
      return e.getStructEntryByStructIndex(structIndex[1:])
    else:
      return e
      
  def getStructEntryByCanonicalIndex(self,indices):
    return self.getStructEntryByStructIndex(filter(lambda x: isinstance(x,str),indices))

  def traverseByPowerIndex(self,powerIndex,canonicalIndex=(),dispatcher=None,payload=None):
      try:
        if len(powerIndex)==0: return dispatcher(payload,canonicalIndex)
        p = powerIndex[0]
        if isinstance(p,str):
          return self.dict[p].traverseByPowerIndex(
            powerIndex[1:],
            canonicalIndex=canonicalIndex+(p,),
            dispatcher=dispatcher,
            payload=payload
          )
        elif isinstance(p,slice):
          raise Exception("slice not allowed here, did you mean '...' ?")
        elif isinstance(p,type(Ellipsis)):
          """
             Why ellipsis? Because it's all or nothing
          """
          return [
                     self.dict[k].traverseByPowerIndex(
                       powerIndex[1:],
                       canonicalIndex=canonicalIndex+(k,),
                       dispatcher=dispatcher,
                       payload=payloadUnpack(payload,i)
                     ) 
                  for i,k in enumerate(self.keys())]
        elif isinstance(p,dict) or isinstance(p,NestedDictLiteral):
          if isinstance(payload,dict):
            return dict([
                    ( k,
                      self.dict[k].traverseByPowerIndex(
                        powerIndex[1:] if isinstance(p,dict) else [p],
                        canonicalIndex=canonicalIndex+(k,),
                        dispatcher=dispatcher,
                        payload=v
                      )
                    ) for k,v in payload.iteritems()
                   ])
          else:
            return dict([
                    ( k,
                      v.traverseByPowerIndex(
                        powerIndex[1:] if isinstance(p,dict) else [p],
                        canonicalIndex=canonicalIndex+(k,),
                        dispatcher=dispatcher,
                        payload=payload
                      )
                    ) for k,v in self.dict.iteritems()
                   ])
        elif isinstance(p,list):
          return [
                     self.traverseByPowerIndex(
                       powerIndex[1:],
                       canonicalIndex=canonicalIndex+(s,),
                       dispatcher=dispatcher,
                       payload=payloadUnpack(payload,i)
                     ) 
                 for i,s in enumerate(p)]
        elif callable(p):
          r = p(self.traverseByPowerIndex(
                powerIndex[1:],
                canonicalIndex=canonicalIndex,
                dispatcher=dispatcher.callableInner(),
                payload=payload
              ))
          return dispatcher.callableOuter(payload,canonicalIndex,extraIndex=None,entry=None,inner=r)
        else:
          raise Exception("I don't know what to do with this: %s" % str(p))
      except Exception as e:
        raise Exception("Error occured in struct context with powerIndex %s, at canonicalIndex %s:\n%s" % (str(powerIndex),str(canonicalIndex),str(e)))
      
# Casadi-dependant Structure framework
    
class Dispatcher:
  def __init__(self,**args):
    for k,v in args.items():
      setattr(self,k,v)
      
  def callableInner(self):
    return self
    
  def callableOuter(self,payload,canonicalIndex,extraIndex=None,entry=None,inner=None):
    return inner
    
#Mixins
class CasadiStructureDerivable:

  def __call__(self,arg=0):
    if isinstance(arg,DMatrix):
      a = arg
    else:
      try:
        a = DMatrix(arg)
      except:
        raise Exception("Call to Structure has weird argument: expecting DMatrix-like")
    if a.shape[0] == 1 and a.shape[1] == 1 and self.size>1:
      a = DMatrix.ones(self.size,1)*a
    return DMatrixStruct(self,data=a)
    
  def repeated(self,arg=0):
    if isinstance(arg,DMatrix):
      a = arg
    else:
      try:
        a = DMatrix(arg)
      except:
        raise Exception("Call to Structure has weird argument: expecting DMatrix-like")
    if not(a.shape[1] == self.size):
       raise Exception("Expecting n x %d DMatrix. Got %s" % (self.size,a.dimString()))
    s = struct([entry("t",struct=self,repeat=a.shape[0])])
    numbers = DMatrixStruct(s,data=a,dataVectorCheck=False)
    return numbers.prefix["t"]
    
  def squared(self,arg=0):
    if isinstance(arg,DMatrix):
      a = arg
    else:
      try:
        a = DMatrix(arg)
      except:
        raise Exception("Call to Structure has weird argument: expecting DMatrix-like")
    if a.shape[0] == 1 and a.shape[1] == 1 and self.size>1:
       a = DMatrix.ones(self.size,self.size)*a
    if not(a.shape[1] == a.shape[0] and a.shape[0]==self.size):
       raise Exception("Expecting square DMatrix of size %s. Got %s" % (self.size,a.dimString()))
    s = struct([entry("t",shapestruct=(self,self))])
    numbers = DMatrixStruct(s,data=a,dataVectorCheck=False)
    return numbers.prefix["t"]

  def squared_repeated(self,arg=0):
    if isinstance(arg,DMatrix):
      a = arg
    else:
      try:
        a = DMatrix(arg)
      except:
        raise Exception("Call to Structure has weird argument: expecting DMatrix-like")
    if not(a.shape[1]==self.size and a.shape[0] % self.size == 0):
       raise Exception("Expecting N x square DMatrix. Got %s" % (self.size,a.dimString()))
    s = struct([entry("t",shapestruct=(self,self),repeat=a.shape[0] / self.size)])
    numbers = DMatrixStruct(s,data=a,dataVectorCheck=False)
    return numbers.prefix["t"]
    
class GetterDispatcher(Dispatcher):
  def __call__(self,payload,canonicalIndex,extraIndex=None,entry=None):
    if canonicalIndex in self.struct.map:
      i = performExtraIndex(self.struct.map[canonicalIndex],extraIndex=extraIndex,entry=entry)
      try:
        return self.master[i]
      except Exception as e:
        raise Exception("Error in powerIndex slicing for canonicalIndex %s: %s" % (str(canonicalIndex),str(e)))
    else:
      raise Exception("Canonical index %s does not exist." % str(canonicalIndex))

class SetterDispatcher(Dispatcher):
  def __call__(self,payload,canonicalIndex,extraIndex=None,entry=None):
    if canonicalIndex in self.struct.map:
      i = performExtraIndex(self.struct.map[canonicalIndex],extraIndex=extraIndex,entry=entry)
      try:
        self.master[i] = payload
      except NotImplementedError:
        raise CompatibilityException("Error in canonicalIndex slicing for %s: Incompatible types in a[i]=b with a %s and b %s." % (str(canonicalIndex),str(self.master),str(payload)))
      except Exception as e:
        raise Exception("Error in powerIndex slicing for canonicalIndex %s: %s" % (str(canonicalIndex),str(e)))
    else:
      raise Exception("Canonical index %s does not exist." % str(canonicalIndex))
      
  def callableInner(self):
    return CasadiStructure.IMatrixDispatcher(struct=self.struct)
  
  def callableOuter(self,payload,canonicalIndex,extraIndex=None,entry=None,inner=None):
    try:
      self.master[inner] = payload
    except NotImplementedError:
      raise CompatibilityException("Error in canonicalIndex slicing for %s: Incompatible types in a[i]=b with a %s and b %s." % (str(canonicalIndex),str(self.master),str(payload)))
    except Exception as e:
      raise Exception("Error in powerIndex slicing for canonicalIndex %s: %s" % (str(canonicalIndex),str(e)))
      
class MasterGettable:
  @properGetitem
  def __getitem__(self,powerIndex):
    return self.struct.traverseByPowerIndex(powerIndex,dispatcher=GetterDispatcher(struct=self.struct,master=self.master))

class MasterSettable:
  @properGetitem
  def __setitem__(self,powerIndex,value):
    return self.struct.traverseByPowerIndex(powerIndex,dispatcher=SetterDispatcher(struct=self.struct,master=self.master),payload=value)
    
def delegation(extraIndex,entry,i):
  if isinstance(extraIndex,str) or (isinstance(extraIndex,list) and len(extraIndex)>0 and all([isinstance(e,str) for e in extraIndex])):
    extraIndex = FlatIndexDelegater(extraIndex)
  if isinstance(extraIndex,Delegater):
    if entry is None: raise Exception("Cannot use delayed index without supplied entry.")
    if entry.shapestruct is None: raise Exception("Cannot use delayed index without supplied shapestruct.")
    if not(isinstance(entry.shapestruct[i],Structure)) : raise Exception("Cannot use delayed index with a integer shapestruct argument.")
    return extraIndex(entry.shapestruct[i])
  else:
    return extraIndex
    
def performExtraIndex(i,extraIndex=None,entry=None):
  if extraIndex is not None and not(isinstance(extraIndex[0],NestedDictLiteral)):
    if len(extraIndex)>2 or len(extraIndex)==0:
      raise Exception("Powerindex exhausted. Remaining %s is interpreted as extraIndex, but length must be 1 or 2." % str(extraIndex))
    try:
      if len(extraIndex)==1:
        a = extraIndex[0]
        a = delegation(a,entry,0)
        return i.__getitem__(a)
      else: 
        a,b = extraIndex
        a = delegation(a,entry,0)
        b = delegation(b,entry,1)
        return i.__getitem__((a,b))
    except NotImplementedError:
       raise Exception("Powerindex exhausted. Passing on %s to %s, but it doesn't know what to do with it" % (str(extraIndex),str(type(i))))
  else:
    return i


class Prefixer:
  def __init__(self,struct,prefix):
    self.struct = struct
    self.prefix = prefix
    
    methods = [ "__DMatrix__", "__SXMatrix__","__MX__"]
    for m in methods:
      if hasattr(self.struct,m):
        setattr(self,m,self.__call__)
          
 
  def __getattr__(self,name):
    # When attributes are not found, delegate to self()
    # This allows for e.g. sin(x) and x+1 to work
    t = self()
    if not(isinstance(t,list) or isinstance(t,dict) or isinstance(t,tuple)):
      return getattr(t,name)
        
  def __str__(self):
    return "prefix( " + str(self.prefix) + "," + self.struct.__str__(compact=True) + ")"
    
  __repr__ = __str__
  
  def __call__(self):
    return self.struct.__getitem__(self.prefix)
  
  @properGetitem
  def __getitem__(self,powerIndex):
    return self.struct.__getitem__(self.prefix + powerIndex)
    
  @properGetitem
  def __setitem__(self,powerIndex,data):
    return self.struct.__setitem__(self.prefix + powerIndex,data)
    
class PrefixConstructor:

  def __str__(self):
    return "prefixConstructor(" + self.struct.__str__(compact=True) + ")"
    
  __repr__ = __str__
  
  def __init__(self,struct):
    self.struct = struct
  
  @properGetitem
  def __getitem__(self,prefix):
    return Prefixer(self.struct,prefix)
    
class CasadiStructure(Structure,CasadiStructureDerivable):
  """
    size
    map
  """

  class FlatIndexDispatcher(Dispatcher):
    def __call__(self,payload,canonicalIndex,extraIndex=None,entry=None):
      if canonicalIndex in self.struct.map:
        return list(performExtraIndex(self.struct.map[canonicalIndex],extraIndex=extraIndex,entry=entry))
      else:
        raise Exception("Canonical index %s not found." % str(canonicalIndex))

  class IMatrixDispatcher(Dispatcher):
    def __call__(self,payload,canonicalIndex,extraIndex=None,entry=None):
      if canonicalIndex in self.struct.map:
        return performExtraIndex(self.struct.map[canonicalIndex],extraIndex=extraIndex,entry=entry)
      else:
        raise Exception("Canonical index %s not found." % str(canonicalIndex))
        
  def __init__(self,*args,**kwargs):
    Structure.__init__(self,*args,**kwargs)
    
    self.map = {}
    self.lookuptable = []
    
    hmap = {}
    k = 0 # Global index counter
    for i in self.traverseCanonicalIndex():
      e = self.getStructEntryByCanonicalIndex(i)
      sp = sp_dense(1,1) if e.sparsity is None else e.sparsity
      m = IMatrix(sp,range(k,k+sp.size()))
      k += sp.size()
      it = tuple(i)
      self.map[it] = m
      self.lookuptable+=[(it,kk,p) for kk,p in enumerate(zip(sp.col(),sp.getRow()))]
      for a in canonicalIndexAncestors(it)[1:]:
        if a in hmap:
          hmap[a].append(m)
        else:
          hmap[a] = [m]
    self.size = k
    for k,v in hmap.iteritems():
      hmap[k] = vecNZcat(v)
    
    self.map.update(hmap)
    
    class StructureGetter:
      def __init__(self,struct):
        self.struct = struct
    
    class IMatrixGetter(StructureGetter):
      @properGetitem
      def __getitem__(self,powerIndex):           
        return self.struct.traverseByPowerIndex(powerIndex,dispatcher=CasadiStructure.IMatrixDispatcher(struct=self.struct))

    class FlatIndexGetter(StructureGetter):
      @properGetitem
      def __getitem__(self,powerIndex):
        return flatten(self.struct.traverseByPowerIndex(powerIndex,dispatcher=CasadiStructure.FlatIndexDispatcher(struct=self.struct)))
            
    self.i = IMatrixGetter(self)
    self.f = FlatIndexGetter(self)
    self.struct = self

  def __str__(self,compact=False):
    return ("" if compact else "Structure with total size %d.\n" % self.size)+ Structure.__str__(self,compact=compact)
    
  def getCanonicalIndex(self,i,extraMode=1):
    """
      Returns the canonicalIndex of the entry with a given flatIndex
      extraMode influences wether nothing (0), [i] (1) or [i,j] (2) will be returned as extra index
    """
    if i<0 or i>=self.size:
      raise Exception("Lookup index out of range. Got %d, but structure is of size %d" % (i,self.size)) 
    can,k,p = self.lookuptable[i]     
    if extraMode==0:
      return can
    elif extraMode==1:
      return can+(k,)
    else:
      return can+p
      
  def canonicalIndices(self,extraMode=1):
    return [self.getCanonicalIndex(i,extraMode=extraMode) for i in range(self.size)]
      
  def getLabel(self,i,extraMode=1):
    t = self.getCanonicalIndex(i,extraMode=extraMode)
    return "["+ ",".join(map(str,t)) + "]"
    
  def labels(self,extraMode=1):
    return [self.getLabel(i,extraMode=extraMode) for i in range(self.size)]
    
class Structured:
  description = "Generic Structured object"
  def __init__(self,structure):
    self.struct = structure.struct
    self.i = self.struct.i
    self.f = self.struct.f
    self.prefix = PrefixConstructor(self)
       
  @property
  def size(self):
    return self.struct.size
    
  @property
  def cat(self):
    return self.master
    
  def __str__(self,compact=False):
    if compact is False:
      return self.description + " with following structure:\n" + self.struct.__str__()
    else:
      return self.description + " (" + self.struct.__str__(compact=True) + ")"
    
  def keys(self):
    return self.struct.keys()
    
class CasadiStructured(Structured,CasadiStructureDerivable):
  description = "Generic Structured object"
  
  def __init__(self,struct,order=None):
    if hasattr(struct,"struct"):
      Structured.__init__(self,struct.struct)
      self.entries = []
    else:
      entrylist = EntryList(struct,order=order)
      self.entries = entrylist.entries
      Structured.__init__(self,CasadiStructure(self.entries, order=entrylist.order))
  
    self.getCanonicalIndex = self.struct.getCanonicalIndex
    self.canonicalIndices = self.struct.canonicalIndices
    self.getLabel = self.struct.getLabel
    self.labels = self.struct.labels
      
  @property
  def shape(self):
    return (self.size,1)

  def sparsity(self):
    return sp_dense(self.size,1)
    
  def getCanonicalIndex(self,*args,**kwargs):
    return self.struct.lookup(*args,**kwargs)

class CompatibilityException(Exception):
  pass

class ssymStruct(CasadiStructured,MasterGettable):
  description = "ssym"
  def __init__(self,struct,order=None):
    CasadiStructured.__init__(self,struct,order=order)
    
    if any(e.expr is not None for e in self.entries):
      raise Exception("struct_ssym does not accept entries with an 'expr' argument, because such an element is not purely symbolic.")
      
    s = []
    for i in self.struct.traverseCanonicalIndex():
      e = self.struct.getStructEntryByCanonicalIndex(i)
      s.append(ssym("_".join(map(str,i)),e.sparsity.size()))
        
    self.master = vecNZcat(s)

    for e in self.entries:
      if e.sym is not None:
        self.master[self.i[e.name]] = e.sym

  def __SXMatrix__(self):
    return self.cat
    
class msymStruct(CasadiStructured,MasterGettable):
  description = "msym"
  def __init__(self,struct,order=None):
    CasadiStructured.__init__(self,struct,order=order)

    if any(e.expr is not None for e in self.entries):
      raise Exception("struct_msym does not accept entries with an 'expr' argument, because such an element is not purely symbolic.")
    if any(e.sym is not None for e in self.entries):
      raise Exception("struct_msym does not accept entries with an 'sym' argument.")

    self.master = msym("V",self.size,1)

  def __MX__(self):
    return self.cat
    
    


class MatrixStruct(CasadiStructured,MasterGettable,MasterSettable):

  @property
  def description(self):
    return "Mutable " + self.mtype.__name__
    
  def __init__(self,struct,mtype,data=None,order=None,dataVectorCheck=True):
    CasadiStructured.__init__(self,struct,order=None)
    
    if any(e.expr is None for e in self.entries):
      raise Exception("struct_SX does only accept entries with an 'expr' argument.")

    self.mtype = mtype
    if isinstance(data,mtype):
      self.master = data
    elif data is None:
      self.master = mtype.nan(self.size,1)
    else:
      self.master = mtype(data)
      
    if dataVectorCheck:
      if self.master.shape[0]!=self.size:
        raise Exception("MatrixStruct: dimension error. Expecting %d-by-1, but got %s" % (self.size,self.master.dimString()))
      if self.master.shape[1]!=1:
        raise Exception("MatrixStruct: dimension error. Expecting %d-by-1, but got %s" % (self.size,self.master.dimString()))
    else:
      if self.master.size()!=self.size:
        raise Exception("MatrixStruct: dimension error. Expecting %d entries, but got %s" % (self.size,self.master.dimString()))
    for e in self.entries:
      self[e.name] = e.expr
      
class DMatrixStruct(MatrixStruct):
  def __init__(self,struct,data=None,dataVectorCheck=True):
    MatrixStruct.__init__(self,struct,DMatrix,data=data,dataVectorCheck=dataVectorCheck)
    
  def __DMatrix__(self):
    return self.cat

class SXMatrixStruct(MatrixStruct):
  def __init__(self,struct,data=None):
    MatrixStruct.__init__(self,struct,SXMatrix,data=data)

  def __SXMatrix__(self):
    return self.cat
    
class MXStruct(MatrixStruct):
  def __init__(self,struct,data=None):
    MatrixStruct.__init__(self,struct,MX,data=data)

  def __MX__(self):
    return self.cat

class MXVeccatStruct(CasadiStructured,MasterGettable):
  description = "Partially mutable MX"
  def __init__(self,arg,order=None):
    CasadiStructured.__init__(self,arg,order=order)
    if any(e.expr is None for e in self.entries):
      raise Exception("struct_MX does only accept entries with an 'expr' argument.")

    self.storage = []
    self.mapping = {}
    for k,i in enumerate(self.struct.traverseCanonicalIndex(limit=1)):
      self.storage.append(None)
      self.mapping[tuple(i)] = k
      
    for e in self.entries:
      self[e.name] = e.expr
      
    self.dirty = True

  def __setitem__(self,powerIndex,value):
    if not isinstance(powerIndex,tuple):
      powerIndex = (powerIndex,)    
        
    def inject(payload,canonicalIndex,extraIndex=None,entry=None):
      if extraIndex is not None:
        raise Exception("An MX veccat structure does not accept indexing on MX level for __setitem__.")
      if not hasattr(self,"sparsity"):
        raise Exception("An MX veccat structure __setitem__ accepts only objects that have sparsity.")
      
      if canonicalIndex in self.mapping:
        if self.struct.map[canonicalIndex].sparsity()!=payload.sparsity():
          raise Exception("Error in powerIndex slicing %s for canonicalIndex %s: Shape mismatch. lhs is %s, rhs is %s." % (str(powerIndex),str(canonicalIndex),self.struct.map[canonicalIndex].sparsity().dimString(),payload.sparsity().dimString()))
        self.storage[self.mapping[canonicalIndex]] = payload
      else:
        raise Exception("Not found: %s " % str(canonicalIndex))
      self.dirty = True
    return self.struct.traverseByPowerIndex(powerIndex,dispatcher=inject,payload=value)
    
  def __MX__(self):
    return self.cat

  @property
  def master(self):
    if any(e is None for e in self.storage):
      missing = filter(lambda k: self.storage[self.mapping[k]] is None,self.mapping.keys())
      
      raise Exception("Problem in MX vecNZcat structure cat: missing expressions. The following entries are missing: %s" % str(missing))
      
    if self.dirty:
      self.master_cached = vecNZcat(self.storage)

    return self.master_cached
    
    
struct_ssym = ssymStruct
struct_msym = msymStruct
struct_SX = SXMatrixStruct
struct_MX_mutable = MXStruct
struct_MX = MXVeccatStruct
struct = CasadiStructured



entry = StructEntry

class CasadiStructEntry(StructEntry):
  def __init__(self,*args,**kwargs):
    if len(args)==0:
      raise Exception("Missing name argument (first argument of Entry)")
    else:
      self.name = args[0]
    self.dict = kwargs
    
    if len(args)>1:
      raise Exception("Don't know what to do with unnamed arguments %s" % str(args[1:]))
      
    

    kw = kwargs.keys()
    kws = ['repeat','shape','sym','expr','struct','shapestruct']
    for k in kw:
      if k not in kws:
        raise Exception("Unknown keyword argument '%s'. Please use one of %s." % (k,str(kws)))
    
    
    
    for kc, fk in [
          ('shape',['struct']),
          ('struct',['shape','shapestruct']),
          ('shapestruct',['struct']), # You might have a sparse matrix with shapestruct
          ('sym',['shape','repeat','expr']),
          ('expr',['shape','repeat','sym'])
        ]:
        if kc in kwargs:
          for fki in fk:
            if fki in kwargs:
              raise Exception("You supplied keyword argument '%s', but it cannot be combined with keyword argument '%s'." % (kc,fki))
    
    #     repeat   argument
    self.repeat = []
    
    if 'repeat' in kwargs:
      self.repeat = kwargs["repeat"] if isinstance(kwargs["repeat"],list) else [kwargs["repeat"]]
    
    if not all(map(lambda x: isinstance(x,int),self.repeat)):
      raise Exception("The 'repeat' argument, if present, must be a list of integers, but got %s" % str(self.repeat))

      
    self.struct = None
    #     struct   argument
    if 'struct' in kwargs:
      struct = kwargs["struct"]
      if isinstance(struct,Structure):
        self.struct = struct
      elif isinstance(struct,Structured):
        self.struct = struct.struct
        
    
    self.sparsity = None
    #     shape   argument
    if 'shape' in kwargs:
      shape = kwargs["shape"]
      if isinstance(shape,int):
        self.sparsity = sp_dense(shape,1)
      elif isinstance(shape,list) or isinstance(shape,tuple):
        if len(shape)==0 or len(shape)>2:
          raise Exception("The 'shape' argument, if present, must be an integer, a tuple of 1 or 2 integers, a sparsity pattern.")
        else:
          self.sparsity = sp_dense(*shape)
      elif isinstance(shape,CRSSparsity):
        self.sparsity = shape
      else:
        raise Exception("The 'shape' argument, if present, must be an integer, a tuple of 1 or 2 integers, or a sparsity pattern. Got %s " % str(shape))
    else:
      self.sparsity = sp_dense(1,1)
    
    self.shapestruct = None
    #     shapestruct  argument
    if 'shapestruct' in kwargs:
      shapestruct = kwargs["shapestruct"]
      if isinstance(shapestruct,Structured) or isinstance(shapestruct,Structure):
        self.shapestruct = (shapestruct.struct,1)
      elif isinstance(shapestruct,tuple):
        if not(all([isinstance(e,Structured) or isinstance(e,Structure) or isinstance(e,int) for e in shapestruct])) or len(shapestruct)==0 or len(shapestruct)>2:
          raise Exception("The 'shapestruct' argument, if present, must be a structure or a tuple of structures or numbers")
        self.shapestruct = tuple([e if isinstance(e,int) else e.struct for e in shapestruct])
      else:
        raise Exception("The 'shapestruct' argument, if present, must be a structure or a tuple of at most structures")
      
      if 'shape' not in kwargs:
        self.sparsity = sp_dense(*[e if isinstance(e,int) else e.size for e in self.shapestruct])
        
    #     sym    argument
    self.sym = None
    if 'sym' in kwargs:
      sym = kwargs["sym"]
      if isinstance(sym,SXMatrix) and isSymbolicSparse(sym):
        self.sym = sym
      elif isinstance(sym,Structured): 
        self.struct = sym.struct
        self.sym = sym.cat
      else:
        raise Exception("The 'sym' argument must be a purely symbolic SXMatrix or a structured symbolic. Got %s instead." % str(self.sym)) 
      self.sparsity = self.sym.sparsity()
      
    #     expr    argument
    self.expr = None
    if 'expr' in kwargs:
      self.expr = kwargs["expr"]
      
      def getPrimitive(e,repeat=[]):
        if isinstance(e,list):
          if len(e)==0:
            return None,repeat+[0]
          else:
            return getPrimitive(e[0],repeat=repeat+[len(e)])
        else:
          return e,repeat
      
          
      p,r = getPrimitive(self.expr)
      
      self.repeat = r
      
      if hasattr(p,"sparsity"):
        self.sparsity = p.sparsity()
      else:
        raise Exception("The 'expr' argument must be a matrix expression or nested list of matrix expressions. Got %s instead." % str(p)) 
      
    StructEntry.__init__(self,self.name,struct=self.struct,dims=self.repeat,data=self.sparsity)
 
  def primitiveString(self):
    return self.sparsity.dimString()
 
def entry(*args,**kwargs):
  if len(args)==1 and isinstance(args[0],CasadiStructEntry):
    return args[0]
  return CasadiStructEntry(*args,**kwargs) 

class EntryList:
  def __init__(self,arg,order = None):
    self.entries = []
    self.order = []
    
    if not isinstance(arg,list):
      raise Exception("Expecting list of entries, with possible tuples for grouping, but got %s" % str(arg))
    
    for e in arg:
      if isinstance(e,tuple):
        entries = map(entry,e)
        self.order.append(tuple(x.name for x in entries))
        self.entries+=entries
      else:
        ee = entry(e)
        self.order.append(ee.name)
        self.entries.append(ee)

    # Override order
    if order is not None:
      if any(isinstance(e,tuple) for e in self.order):
        raise Exception("You supplied an order by using tuple syntax on entries %s, but you overwrite it with the 'order' keyword. Use one or the other, not both.")
      self.order = order
      
    self.names = map(lambda x : x.name,self.entries)
    if len(self.names)!=len(set(self.names)):
      duplicates = []
      for i,e in enumerate(self.names):
        if e in self.names[:i] or e in self.names[i+1:]:
          duplicates.append(e)
      raise Exception("Your list of entries contains duplicates: %s" % str(list(set(duplicates))))


class Delegater:
  def __init__(self,arg):
    self.arg = arg

  def __str__(self):
    return "%s[%s]" % (self.__class__.__name__,str(self.arg))
    
  __repr__ = __str__


class IndexDelegater(Delegater):
  def __call__(self,struct):
    return struct.i.__getitem__(self.arg)

class FlatIndexDelegater(Delegater):
  def __call__(self,struct):
    return struct.f.__getitem__(self.arg)

    
class DelegaterConstructor:
  """
    Creates an object that delegates a slicing operation.
    
    Example usage:
      s = struct_ssym([])
      x = struct_ssym(entry("x",sp_diag(4)))
      x["x",0,index[:]]
    
  """
  def __init__(self,delegater,prepend=()):
    self.prepend = prepend
    self.delegater = delegater
    
  @properGetitem
  def __getitem__(self,arg):
    return self.delegater(self.prepend + arg)
    
index  = DelegaterConstructor(IndexDelegater)
indexf = DelegaterConstructor(FlatIndexDelegater)
