from casadi import *

import operator

# StructIndex :tuple/list of strings
# canonicalIndex : tuple/list of string or numbers
# powerIndex: tuple/list of string, numbers, lists, slices, dicts

def lpack(L): return [[x] for x in L]
  
def combine(*args):
  if len(args)==0:
    return [[]]
  else:
    return [a + b for a in args[0] for b in combine(*args[1:])]

def listindices(dims,nest=False):
  if len(dims)==0:
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

class Repeater:
  def __init__(self,e): self.e = e

def repeat(e):  return Repeater(e)

def canonical(ind,s):
  if ind < 0:
    return ind + s
  else:
    return ind

class StructEntry:
  def __init__(self,name,struct=None,data=None,dims=[]):
    self.name = name
    self.data = data
    self.dims = dims
    self.struct = struct

  def __str__(self,compact=False):
     s=''
     if len(self.dims)==1:
       s+= "list(%d): " % self.dims[0]
     if len(self.dims)>1:
       s+= "nestedlist(%s): " % str(self.dims)
     if self.isPrimitive():
       if hasattr(self.data,'dimString'):
         s+= self.data.dimString()
       else:
         s+= "*"
     else:
       s+=self.struct.__str__(compact=True)
     return s
     
  def isPrimitive(self):
    return self.struct is None
    
  def traverseCanonicalIndex(self,nest=False):
    children = [[]] if self.struct is None else self.struct.traverseCanonicalIndex()
    li = listindices(self.dims,nest)
    n = [[self.name]]
    if nest:
      return [combine(n,i,children) for i in li]
    else:
      return combine(n,li,children)
      
  def getStructEntryByStructIndex(self,structIndex):
    return self.struct.getStructEntryByStructIndex(structIndex)
      
  def traverseByPowerIndex(self,powerIndex,dims=None,canonicalIndex=(),dispatcher=None,payload=None):
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
        return dispatcher(payload,canonicalIndex)
        
    if len(dims)==0:
      if self.isPrimitive(): # Pass on remainder of powerIndex to dispatcher
        return dispatcher(payload,canonicalIndex,tuple(powerIndex))
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
        if isinstance(payload,list):
          newpayloads = payload
        elif isinstance(payload,Repeater):
          newpayloads = [payload.e] * s
        else:
          newpayloads = [payload] * s
        return [
                  self.traverseByPowerIndex(
                    powerIndex[1:],
                    dims=dims[1:],
                    canonicalIndex=canonicalIndex+(canonical(i,s),),
                    dispatcher=dispatcher,
                    payload = newpayloads[i]
                  )
               for i in p]
      elif callable(p):
        r = self.traverseByPowerIndex(
              powerIndex[1:],
              dims=dims,
              canonicalIndex=canonicalIndex,
              dispatcher=dispatcher,
              payload=payload
            )
        return p(r)
      else:
        raise Exception("I don't know what to do with this: %s" % str(p))
        
class Structure:
  def __init__(self,entries,order=None):
    self.entries = entries
    self.order = [e.name for e in self.entries] if order is None else order
    self.dict = dict([(e.name,e) for e in self.entries])
    
    for e in self.order:
      if isinstance(e,str):
        if e not in self.dict:
          raise Exception("Order '%s' is invalid." % e)
      elif isinstance(e,tuple):
        for ee in e:
          if ee not in self.dict:
            raise Exception("Order '%s' is invalid." % ee)
       

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
     
     
        
  def traverseCanonicalIndex(self):
    ret = []
    for d in self.order:
      if isinstance(d,tuple):
        for v in intersperse(*[self.dict[de].traverseCanonicalIndex(True) for de in d]):
          ret += v 
      else:
        ret += self.dict[d].traverseCanonicalIndex()
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
      raise Exception("Slice not permitted here")
    elif isinstance(p,dict):
      if isinstance(payload,dict):
        return dict([
                ( k,
                  self.dict[k].traverseByPowerIndex(
                    powerIndex[1:],
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
                    powerIndex[1:],
                    canonicalIndex=canonicalIndex+(k,),
                    dispatcher=dispatcher,
                    payload=payload
                  )
                ) for k,v in self.dict.iteritems()
               ])
    elif isinstance(p,list):
      if isinstance(payload,list):
        return [
                 self.traverseByPowerIndex(
                   powerIndex[1:],
                   canonicalIndex=canonicalIndex+(s,),
                   dispatcher=dispatcher,
                   payload=payload[i]
                 ) 
               for i,s in enumerate(p)]
      elif isinstance(payload,Repeater):
        return [
                 self.traverseByPowerIndex(
                   powerIndex[1:],
                   dims[1:0],
                   canonicalIndex=canonicalIndex+(s,),
                   dispatcher=dispatcher,
                   payload=payload.e
                 )
               for s in p]
      else:
        return [
                 self.traverseByPowerIndex(
                   powerIndex[1:],
                   canonicalIndex=canonicalIndex+(s,),
                   dispatcher=dispatcher,
                   payload=payload
                 )
               for s in p]
    elif callable(p):
      r = self.traverseByPowerIndex(
            powerIndex[1:],
            canonicalIndex=canonicalIndex,
            dispatcher=dispatcher,
            payload=payload
          )
      return p(r)
    else:
      raise Exception("I don't know what to do with this: %s" % str(p))
      
class CasadiStructureDerivable:

  def __call__(self,arg=None):
    if isinstance(arg,DMatrix):
      return self.DMatrix(arg)
    else:
      raise Exception("Call to Structure has weird argument: expecting DMatrix")
    
  def ssym(self):
    return ssymStruct(self)

  def msym(self,createParent=True):
    return msymStruct(self,createParent=createParent)
    
  def zeros(self):
    return DMatrixStruct(self)

  def DMatrix(self,data=None):
    return DMatrixStruct(self,data=data)

  def SXMatrix(self,data=None):
    return SXMatrixStruct(self,data=data)

  def MX(self,data=None):
    return MXStruct(self,data=data)

def performExtraIndex(i,extraIndex):
  if extraIndex is not None:
    if len(extraIndex)==1:
      return i.__getitem__(*extraIndex)
    else: 
      return i.__getitem__(extraIndex)
  else:
    return i
      
class CasadiStructure(Structure,CasadiStructureDerivable):
  """
    size
    map
  """

  def __init__(self,*args,**kwargs):
    Structure.__init__(self,*args,**kwargs)
    
    self.map = {}
    
    hmap = {}
    k = 0 # Global index counter
    for i in self.traverseCanonicalIndex():
      e = self.getStructEntryByCanonicalIndex(i)
      sp = sp_dense(1,1) if e.data is None else e.data
      m = IMatrix(sp,range(k,k+sp.size()))
      k += sp.size()
      it = tuple(i)
      self.map[it] = m
      for a in canonicalIndexAncestors(it)[1:]:
        if a in hmap:
          hmap[a].append(m)
        else:
          hmap[a] = [m]
    self.size = k
    for k,v in hmap.iteritems():
      hmap[k] = veccat(v)
    
    self.map.update(hmap)
    
    class IMatrixGetter:
      def __init__(self,struct):
        self.struct = struct
        
      def __getitem__(self,powerIndex):
        if not isinstance(powerIndex,tuple):
          powerIndex = (powerIndex,)    
            
        def inject(payload,canonicalIndex,extraIndex=None):
          if canonicalIndex in self.struct.map:
            return performExtraIndex(self.struct.map[canonicalIndex],extraIndex)
          
        return self.struct.traverseByPowerIndex(powerIndex,dispatcher=inject)
    self.i = IMatrixGetter(self)

    
class Structured:
  def __init__(self,structure):
    self.struct = struct(structure)
    self.i = self.struct.i
       
  @property
  def size(self):
    return self.struct.size
    
  @property
  def cat(self):
    return self.master
    
class CasadiStructured(Structured,CasadiStructureDerivable):
  @property
  def shape(self):
    return (self.size,1)

class CompatibilityException(Exception):
  pass

class MasterGettable:
  def __getitem__(self,powerIndex):
    if not isinstance(powerIndex,tuple):
      powerIndex = (powerIndex,)    
        
    def inject(payload,canonicalIndex,extraIndex=None):
      if canonicalIndex in self.struct.map:
        i = performExtraIndex(self.struct.map[canonicalIndex],extraIndex)
        return self.master[i]
      
    return self.struct.traverseByPowerIndex(powerIndex,dispatcher=inject)

class MasterSettable:
  def __setitem__(self,powerIndex,value):
    if not isinstance(powerIndex,tuple):
      powerIndex = (powerIndex,)    
        
    def inject(payload,canonicalIndex,extraIndex=None):
      if canonicalIndex in self.struct.map:
        i = performExtraIndex(self.struct.map[canonicalIndex],extraIndex)
        try:
          self.master[i] = payload
        except NotImplementedError:
          raise CompatibilityException("Incompatible types in a[i]=b with a %s and b %s" % (str(type(self.master)),str(type(payload))))
        
    return self.struct.traverseByPowerIndex(powerIndex,dispatcher=inject,payload=value)
    
class ssymStruct(CasadiStructured,MasterGettable):
  def __init__(self,struct):
    Structured.__init__(self,struct)
    s = []
    for i in self.struct.traverseCanonicalIndex():
      e = self.struct.getStructEntryByCanonicalIndex(i)
      sp = sp_dense(1,1) if e.data is None else e.data
      s.append(ssym("_".join(map(str,i)),sp.size()))
    self.master = veccat(s)

  def __SXMatrix__(self):
    return self.cat
    
class msymStruct(CasadiStructured,MasterGettable):
  def __init__(self,struct,createParent=True):
    Structured.__init__(self,struct)
    if createParent:
      self.master = msym("V",self.size,1)
    else:
      s = []
      for i in self.struct.traverseCanonicalIndex():
        e = self.struct.getStructEntryByCanonicalIndex(i)
        sp = sp_dense(1,1) if e.data is None else e.data
        s.append(msym("_".join(map(str,i)),sp.size()))
      self.master = veccat(s)

  def __MX__(self):
    return self.cat
    
class MatrixStruct(CasadiStructured,MasterGettable,MasterSettable):
  def __init__(self,struct,mtype,data=None):
    Structured.__init__(self,struct)
    self.mtype = mtype
    self.master = mtype(data) if data is not None else mtype.zeros(self.size,1)
    if self.master.shape[0]!=self.size:
      raise Exception("MatrixStruct: dimension error. Expecting %d-by-1, but got %s", (self.size,self.master.dimString()))
    if self.master.shape[1]!=1:
      raise Exception("MatrixStruct: dimension error. Expecting %d-by-1, but got %s", (self.size,self.master.dimString()))
      
class DMatrixStruct(MatrixStruct):
  def __init__(self,struct,data=None):
    MatrixStruct.__init__(self,struct,DMatrix,data=data)
    
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
    
states = CasadiStructure([
            StructEntry('x'),
            StructEntry('y'),
            StructEntry('z'),
            StructEntry('u',data=sp_dense(4)),
            StructEntry('v',dims=[4,2]),
            StructEntry('w',dims=[6]),
            StructEntry('p',dims=[9],data=sp_dense(6))
         ],order=['x','y','z','u',('v','w'),'p'])
         
shooting = CasadiStructure([StructEntry('X',struct=states,dims=[4,5]),StructEntry('U',dims=[3])],order=[('X','U')])


         
def struct(arg):
  # Style 1: Struct(['x','y','z'])
  # Style 2: Struct([(4,2,'x',2,3),'y','z'])
  # Style 3: Struct([('y',struct),'z'])
  # Style 4: ordering:  Struct([('x','y'),'z'])
  
  if isinstance(arg,Structure):
    return arg
  elif isinstance(arg,Structured):
    return arg.struct
  
  def parsePrimitive(e):
    if isinstance(e,str):
      return e,[StructEntry(e,data=sp_dense(1,1))]
    elif isinstance(e,tuple):
      l = list(e)
      isstring = lambda x: isinstance(x,str)
      istuple = lambda x: isinstance(x,tuple) and any(map(isstring,x))
      strings = map(isstring,l)
      if sum(strings)>1 or any(map(istuple,l)):
        orders = ()
        entries = []
        for ee in l:
          order,entry = parsePrimitive(ee)
          orders += (order,)
          entries+= entry
        return orders,entries
      else:
        si = strings.index(True)
        sp = None
        struct = None
        if si+1==len(l):
          sp = sp_dense(1,1)
        else:
          if isinstance(l[si+1],CRSSparsity):
            sp = l[si+1]
          elif isinstance(l[si+1],int):
            sp = sp_dense(*l[si+1:])
          elif isinstance(l[si+1],Structure):
            struct = l[si+1]
          elif isinstance(l[si+1],tuple):
            sp = sp_dense(*l[si+1])
          else:
            raise Exception("I do not understand this: %s " % str(l[si+1:]))
        return l[si],[StructEntry(l[si],dims=l[:si],data=sp,struct=struct)]
    else: 
      raise Exception("I do not understand this: %s " % str(e))
  
  order,entries = zip(*map(parsePrimitive,arg))
  entries = sum(entries,[])
  
  return CasadiStructure(entries,order=order)
  
def storage(d,order=None):
  order = sorted(d.keys()) if order is None else order
  s = []
  for k in order:
    v = d[k]
    def createEntry(v,dims=[]):
      if isinstance(v,list):
        return createEntry(v[0],dims=dims+[len(v)])
      else:
        return StructEntry(k,dims=dims,data=v.sparsity())
    s.append(createEntry(v))
  struct = CasadiStructure(s,order=order)
  for method in ['DMatrix','SXMatrix','MX']:
    try:
      b = getattr(struct,method)()
      b[{}] = d
      return b
    except CompatibilityException:
      pass
      
      
if __name__=="__main__":  
  print struct(['x','y','z'])

  print struct(['x',(4,3,'y'),'z'])

  print struct(['x',(5,'y'),('z',9)])

  s = struct(['u','v'])

  print struct(['x',('z',s)])

  print struct([('x','y'),'z'])


  print struct([('x',(4,3,'y')),'z'])


  st = storage({'x': ssym("x")**2,'y': [ssym("a",2)**2,ssym("b",2)**2]})

  print st.cat

  print st['x']
  print st['y']


  #raise Exception("x")
           


      
  d = DMatrix([4,5,6])


  print canonicalIndexAncestors(('a',0,0,0,'b',6,'c',2,3,'d','f',0))

  print lpack([1,2,3])
  print lpack('abc')
  print "foo"
  print combine(lpack('abc'))
  print combine(lpack('abc'),lpack([1,2,3]))
  for i in combine(lpack('abc'),lpack([1,2,3]),lpack([4,5])):
    print i

  print combine([[]],lpack([1,2,3]))
  print combine(lpack('abc'),[[]])

  print "bar"
  print listindices([2,2,3])
  print listindices([2,2,3],True)
  print listindices([2,3,4])

  for e in intersperse([1,2,3],[4,5,6,7],[8,9]):
    print e

  print "traverse"
  for e in states.dict['v'].traverseCanonicalIndex():
    print e
  print "traverse"
  for e in states.dict['v'].traverseCanonicalIndex(True):
    print e


       
  print "abc"     
  for e in states.traverseCanonicalIndex():
    print e
    
  print "shooting"
  for e in shooting.traverseCanonicalIndex():
    print e
    

  print shooting.size

  V = ssymStruct(shooting)

  print V.master
  print V.cat

  print V['X',0,0]

  print V['X',0,0,{}]

  print V['X',0,0,['x','y']]

  print V['X',0,0,'v',vertcat,:,horzcat,:]

  print V['X',0,0,'p',:,::2,:]

  print V['X',0,0,'v',horzcat,:,horzcat,:]
  print V['X',0,0,'v',:,horzcat,:]
  print V['X',0,0,'v',blockcat,:,:]

  print V['X',0,:,'v',blockcat]

  print V['X',0,vertcat,:,'v',blockcat]


  print V['X',:,:,['x','y']] 
  #V['X',:,:,['x','y']] = repeat(repeat([6,5]))

  #raise Exception("")

  print V['X',:,:,'v']
  print V['X',:,:,'v',:]
  print V['X',:,:,'v',:,:]

  print V['X',:,:]

  #raise Exception("")

  #print V['U',:,:,:]


  init = DMatrixStruct(shooting)

  print V['X']
  print V['X',:]
  print V['X',:,:]

  print V['X',:,:,{}] 

  print V['X',:,:,'p'] 
  print V['X',:,:,'p',:] 

  print V['X',0,0,'p',:]



  #  :    vertcat
  #  ...  horzcat
  #  []   leave list
  #  {}   leave dict

  #  vertcat - horzcat  : or any operator taking a list argument
  #  {}

  # But: how to do specific ranges that you veccat?


  # New rules
  # At any point in the traverseByPowerIndex list, you can drop a function.
  # It will be applied to whatever is obtained from the left


  # Right hand sides lists? Will not be interpreted as Matrix.

  # Forbidden 

  # same as
  print V['X',:,:]

  print V['X',vertcat,:,horzcat,:]
  print V['X',blockcat,:,:]


  init['X',0,-1,'p',0] = 2
  print init.cat
  init['X',0,-1,'p',:] = 3
  print init.cat
  #init['X',0,0,'p',:] = range(9)
  #print init.cat()

  init['X',:,:,['x','y']] = repeat(repeat([6,5]))

  print init
  print init.cat

  init['X',:,:,{}] = repeat(repeat({'x': 9,'y': 3}))

  print init
  print init.cat

