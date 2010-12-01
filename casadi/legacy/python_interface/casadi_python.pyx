# -*- coding: utf-8 -*-
cimport casadi_python
cimport python_exc
import numpy as np
cimport numpy as np

ADD_NODE = 0
SUB_NODE = 1
MUL_NODE = 2
DIV_NODE = 3
NEG_NODE = 4
EXP_NODE = 5
LOG_NODE = 6
POW_NODE = 7
SQRT_NODE = 8
SIN_NODE = 9
COS_NODE = 10
TAN_NODE = 11
ASIN_NODE = 12
ACOS_NODE = 13
ATAN_NODE = 14
STEP_NODE = 15
FLOOR_NODE = 16
CEIL_NODE = 17
EQUALITY_NODE = 18
ERF_NODE = 19
FMIN_NODE = 20
FMAX_NODE = 21

cdef class STLString:
  cpdef casadi_python.string_ptr _string_ptr

  def __cinit__(self, *args):
    self._string_ptr = casadi_python.casadi_string_new()
    if self._string_ptr is NULL:
      python_exc.PyErr_NoMemory()
      
    if len(args)>0:
      casadi_string_assign(self._string_ptr, args[0])
      
  def __dealloc__(self):
    if self._string_ptr is not NULL:
      casadi_python.casadi_string_delete(self._string_ptr)

  def __str__(self):
    return casadi_string_get(self._string_ptr)
    
  def __repr__(self):
    return casadi_string_get(self._string_ptr)


cdef class STLVector:
  cpdef casadi_python.vector_ptr _vector_ptr

  def __cinit__(self, *args):
    self._vector_ptr = casadi_python.casadi_vector_new()
    if self._vector_ptr is NULL:
      python_exc.PyErr_NoMemory()
      
  def __dealloc__(self):
    if self._vector_ptr is not NULL:
      casadi_python.casadi_vector_delete(self._vector_ptr)

  cpdef resize(self, int n):
    casadi_vector_resize(self._vector_ptr, n)

  def set_double(self, np.ndarray[np.float64_t, ndim=1] v):
    if self._vector_ptr is not NULL:
      assert casadi_vector_size(self._vector_ptr) == len(v), 'size mismatch'
      casadi_vector_set(self._vector_ptr, <double*>v.data)

  def get_double(self, np.ndarray[np.float64_t, ndim=1] v):
    if self._vector_ptr is not NULL:
      assert casadi_vector_size(self._vector_ptr) == len(v), 'size mismatch'
      casadi_vector_get(self._vector_ptr, <double*>v.data)
    

cdef class MX:
  cpdef casadi_python.mx_ref _mx_ref

  def __cinit__(self, *args):
    self._mx_ref= casadi_python.casadi_mx_new()
    if self._mx_ref is NULL:
      python_exc.PyErr_NoMemory()

    # quick return if default constructor
    if len(args)==0:
      return
    
    # data array
    data = args[0]

    # temp
    cdef double tmp[1]
    cdef char order = ord("R") # row major

    if isinstance(data,str):
      casadi_python.casadi_mx_symbol(self._mx_ref,data,1,1)
    elif isinstance(data,int):
      tmp[0] = data
      casadi_python.casadi_mx_constant(self._mx_ref,tmp,1,1,order)
    elif isinstance(data,float):
      tmp[0] = data
      casadi_python.casadi_mx_constant(self._mx_ref,tmp,1,1,order)
      
  cpdef print_string(self, STLString s):
    casadi_python.casadi_mx_print_string(self._mx_ref,s._string_ptr)

  def __repr__(self):
    s = STLString()
    self.print_string(s)
    return str(s)

  def __str__(self):
    s = STLString() 
    self.print_string(s)
    return str(s)

  cpdef cbinary(self, int op, MX x, MX y):
    casadi_python.casadi_mx_binary(self._mx_ref,op,x._mx_ref,y._mx_ref)

  def binary(self,op,x,y):
    if not isinstance(x,MX):
      x = MX(x)
    if not isinstance(y,MX):
      y = MX(y)
    self.cbinary(op,x,y)

  cpdef cunary(self, int op, MX x):
    casadi_python.casadi_mx_unary(self._mx_ref,op,x._mx_ref)

  def unary(self,op,x):
    if not isinstance(x,MX):
      x = MX(x)
    self.cunary(op,x)

  def __add__(self,other):
    ret = MX()
    ret.binary(ADD_NODE,self,other)
    return ret

  def __sub__(self,other):
    ret = MX()
    ret.binary(SUB_NODE,self,other)
    return ret

  def __mul__(self,other):
    ret = MX()
    ret.binary(MUL_NODE,self,other)
    return ret

  def __truediv__(self,other):
    ret = MX()
    ret.binary(DIV_NODE,self,other)
    return ret

  def __div__(self,other):
    ret = MX()
    ret.binary(DIV_NODE,self,other)
    return ret

  def __neg__(self):
    ret = MX()
    ret.unary(NEG_NODE,self)
    return ret

  def exp(self):
    ret = MX()
    ret.unary(EXP_NODE,self)
    return ret

  def log(self):
    ret = MX()
    ret.unary(LOG_NODE,self)
    return ret

  def __pow__(self,other,arg3):
    ret = MX()
    ret.binary(POW_NODE,self,other)
    return ret

  def sqrt(self):
    ret = MX()
    ret.unary(SQRT_NODE,self)
    return ret

  def sin(self):
    ret = MX()
    ret.unary(SIN_NODE,self)
    return ret

  def cos(self):
    ret = MX()
    ret.unary(COS_NODE,self)
    return ret

  def tan(self):
    ret = MX()
    ret.unary(TAN_NODE,self)
    return ret

  def arcsin(self):
    ret = MX()
    ret.unary(ASIN_NODE,self)
    return ret

  def arccos(self):
    ret = MX()
    ret.unary(ACOS_NODE,self)
    return ret

  def arctan(self):
    ret = MX()
    ret.unary(ATAN_NODE,self)
    return ret

  def step(self):
    ret = MX()
    ret.unary(STEP_NODE,self)
    return ret

  def floor(self):
    ret = MX()
    ret.unary(FLOOR_NODE,self)
    return ret

  def ceil(self):
    ret = MX()
    ret.unary(CEIL_NODE,self)
    return ret

  def erf(self):
    ret = MX()
    ret.unary(ERF_NODE,self)
    return ret

  def fmin(self,other):
    ret = MX()
    ret.binary(FMIN_NODE,self,other)
    return ret

  def fmax(self,other):
    ret = MX()
    ret.binary(FMAX_NODE,self,other)
    return ret

  def __dealloc__(self):
    if self._mx_ref is not NULL:
      casadi_python.casadi_mx_delete(self._mx_ref)




cdef class SX:
  cpdef casadi_python.sx_ptr _sx_ptr

  def __cinit__(self, *args):
    self._sx_ptr= casadi_python.casadi_sx_new()
    if self._sx_ptr is NULL:
      python_exc.PyErr_NoMemory()

    # quick return if default constructor
    if len(args)==0: return
    
    # argument
    data = args[0]
    
    if isinstance(data,str):
      # variable
      casadi_python.casadi_sx_symbol(self._sx_ptr,data)
    else:
      # constant
      val = np.double(data)
      casadi_python.casadi_sx_constant(self._sx_ptr,val)
      
  cpdef print_string(self, STLString s):
    casadi_python.casadi_sx_print(self._sx_ptr,s._string_ptr)

  def __repr__(self):
    s = STLString() 
    self.print_string(s)
    return str(s)

  def __str__(self):
    s = STLString() 
    self.print_string(s)
    return str(s)

  cpdef cbinary(self, int op, SX x, SX y):
    casadi_python.casadi_sx_binary(self._sx_ptr,op,x._sx_ptr,y._sx_ptr)

  def binary(self,op,x,y):
    if not isinstance(x,SX):
      x = SX(x)
    if not isinstance(y,SX):
      y = SX(y)
    self.cbinary(op,x,y)

  cpdef cunary(self, int op, SX x):
    casadi_python.casadi_sx_unary(self._sx_ptr,op,x._sx_ptr)

  def unary(self,op,x):
    if not isinstance(x,SX):
      x = SX(x)
    self.cunary(op,x)

  def __add__(self,other):
    ret = SX()
    ret.binary(ADD_NODE,self,other)
    return ret

  def __sub__(self,other):
    ret = SX()
    ret.binary(SUB_NODE,self,other)
    return ret

  def __mul__(self,other):
    ret = SX()
    ret.binary(MUL_NODE,self,other)
    return ret

  def __truediv__(self,other):
    ret = SX()
    ret.binary(DIV_NODE,self,other)
    return ret

  def __div__(self,other):
    ret = SX()
    ret.binary(DIV_NODE,self,other)
    return ret

  def __neg__(self):
    ret = SX()
    ret.unary(NEG_NODE,self)
    return ret

  def exp(self):
    ret = SX()
    ret.unary(EXP_NODE,self)
    return ret

  def log(self):
    ret = SX()
    ret.unary(LOG_NODE,self)
    return ret

  def __pow__(self,other,arg3):
    ret = SX()
    ret.binary(POW_NODE,self,other)
    return ret

  def sqrt(self):
    ret = SX()
    ret.unary(SQRT_NODE,self)
    return ret

  def sin(self):
    ret = SX()
    ret.unary(SIN_NODE,self)
    return ret

  def cos(self):
    ret = SX()
    ret.unary(COS_NODE,self)
    return ret

  def tan(self):
    ret = SX()
    ret.unary(TAN_NODE,self)
    return ret

  def arcsin(self):
    ret = SX()
    ret.unary(ASIN_NODE,self)
    return ret

  def arccos(self):
    ret = SX()
    ret.unary(ACOS_NODE,self)
    return ret

  def arctan(self):
    ret = SX()
    ret.unary(ATAN_NODE,self)
    return ret

  def step(self):
    ret = SX()
    ret.unary(STEP_NODE,self)
    return ret

  def floor(self):
    ret = SX()
    ret.unary(FLOOR_NODE,self)
    return ret

  def ceil(self):
    ret = SX()
    ret.unary(CEIL_NODE,self)
    return ret

  def erf(self):
    ret = SX()
    ret.unary(ERF_NODE,self)
    return ret

  def fmin(self,other):
    ret = SX()
    ret.binary(FMIN_NODE,self,other)
    return ret

  def fmax(self,other):
    ret = SX()
    ret.binary(FMAX_NODE,self,other)
    return ret

  def __dealloc__(self):
    if self._sx_ptr is not NULL:
      casadi_python.casadi_sx_delete(self._sx_ptr)


cdef class SXVector:
  cpdef casadi_python.sx_vec_ptr _sx_vec_ptr

  def __cinit__(self, *args):
    self._sx_vec_ptr = casadi_python.casadi_sx_vec_new()
    if self._sx_vec_ptr is NULL:
      python_exc.PyErr_NoMemory()
      
    # quick return if default constructor  
    if len(args)==0: return
    assert len(args)==1
    arg = args[0]

    if isinstance(arg,list) or isinstance(arg,tuple):
      for x in arg:
        if isinstance(x,SX):
          self.push_back(x)
        else:
          self.push_back(SX(x))
    else:
      if isinstance(x,SX):
        self.push_back(x)
      else:
        self.push_back(SX(x))
        
      print "self = ", self
    
  def __dealloc__(self):
    if self._sx_vec_ptr is not NULL:
      casadi_python.casadi_sx_vec_delete(self._sx_vec_ptr)

  cpdef push_back(self, SX ref):
    casadi_sx_vec_push_back(self._sx_vec_ptr, ref._sx_ptr)
    
  cpdef print_string(self, STLString s):
    casadi_python.casadi_sx_vec_print(self._sx_vec_ptr,s._string_ptr)

  def __repr__(self):
    s = STLString()
    self.print_string(s)
    return str(s)

  def __str__(self):
    s = STLString() 
    self.print_string(s)
    return str(s)

















cdef class SXMatrix:
  cpdef casadi_python.sx_matrix_ref _sx_matrix_ref


  def __cinit__(self, *args):
    self._sx_matrix_ref= casadi_python.casadi_sx_matrix_new()
    if self._sx_matrix_ref is NULL:
      python_exc.PyErr_NoMemory()

    # quick return if default constructor
    if len(args)==0:
      return
    
    # data array
    data = args[0]

    # Get dimensions and cow/col major
    cdef int n = 1
    cdef int m = 1
    cdef char order = ord("R") # row major

    if len(args)>=2:
      n = args[1]
    if len(args)>=3:
      n = args[2]
    if len(args)>=4:
      order = ord(args[3])
    assert len(args) <= 4
    
    # temp
    cdef double tmp[1]

    if isinstance(data,str):
      casadi_python.casadi_sx_matrix_symbol(self._sx_matrix_ref,data,n,m)
    elif isinstance(data,int):
      tmp[0] = data
      casadi_python.casadi_sx_matrix_constant(self._sx_matrix_ref,tmp,n,m,order)
    elif isinstance(data,float):
      tmp[0] = data
      casadi_python.casadi_sx_matrix_constant(self._sx_matrix_ref,tmp,n,m,order)
    elif isinstance(data,SX):
      self.sx_matrix_sx(data)
    elif isinstance(data,np.ndarray):
      assert data.ndim==1 or data.ndim==2
      if data.ndim==2: raise Exception('not implemented')
      v = SXVector()
      for x in data:
        if isinstance(x,SX):
          v.push_back(x)
        else:
          v.push_back(SX(x))
      self.sx_matrix_sx_vec(v)
          
    else:
      raise Exception('Cannot convert to SXMatrix',data)
      
  cpdef sx_matrix_sx(self, SX scalar):
    casadi_python.casadi_sx_matrix_sx(self._sx_matrix_ref,scalar._sx_ptr)

  cpdef sx_matrix_sx_vec(self, SXVector v):
    casadi_python.casadi_sx_matrix_sx_vec(self._sx_matrix_ref,v._sx_vec_ptr)

  cpdef print_string(self, STLString s):
    casadi_python.casadi_sx_matrix_print(self._sx_matrix_ref,s._string_ptr)

  def __repr__(self):
    s = STLString() 
    self.print_string(s)
    return str(s)

  def __str__(self):
    s = STLString() 
    self.print_string(s)
    return str(s)

  cpdef cbinary(self, int op, SXMatrix x, SXMatrix y):
    cdef int flag
    flag = casadi_python.casadi_sx_matrix_binary(self._sx_matrix_ref,op,x._sx_matrix_ref,y._sx_matrix_ref)
    assert flag==0

  def binary(self,op,x,y):
    if not isinstance(x,SXMatrix):
      x = SXMatrix(x)
    if not isinstance(y,SXMatrix):
      y = SXMatrix(y)
    self.cbinary(op,x,y)

  cpdef cunary(self, int op, SXMatrix x):
    casadi_python.casadi_sx_matrix_unary(self._sx_matrix_ref,op,x._sx_matrix_ref)

  def unary(self,op,x):
    if not isinstance(x,SXMatrix):
      x = SXMatrix(x)
    self.cunary(op,x)

  def __add__(self,other):
    ret = SXMatrix()
    ret.binary(ADD_NODE,self,other)
    return ret

  def __sub__(self,other):
    ret = SXMatrix()
    ret.binary(SUB_NODE,self,other)
    return ret

  def __mul__(self,other):
    ret = SXMatrix()
    ret.binary(MUL_NODE,self,other)
    return ret

  def __truediv__(self,other):
    ret = SXMatrix()
    ret.binary(DIV_NODE,self,other)
    return ret

  def __div__(self,other):
    ret = SXMatrix()
    ret.binary(DIV_NODE,self,other)
    return ret

  def __neg__(self):
    ret = SXMatrix()
    ret.unary(NEG_NODE,self)
    return ret

  def exp(self):
    ret = SXMatrix()
    ret.unary(EXP_NODE,self)
    return ret

  def log(self):
    ret = SXMatrix()
    ret.unary(LOG_NODE,self)
    return ret

  def __pow__(self,other,arg3):
    ret = SXMatrix()
    ret.binary(POW_NODE,self,other)
    return ret

  def sqrt(self):
    ret = SXMatrix()
    ret.unary(SQRT_NODE,self)
    return ret

  def sin(self):
    ret = SXMatrix()
    ret.unary(SIN_NODE,self)
    return ret

  def cos(self):
    ret = SXMatrix()
    ret.unary(COS_NODE,self)
    return ret

  def tan(self):
    ret = SXMatrix()
    ret.unary(TAN_NODE,self)
    return ret

  def arcsin(self):
    ret = SXMatrix()
    ret.unary(ASIN_NODE,self)
    return ret

  def arccos(self):
    ret = SXMatrix()
    ret.unary(ACOS_NODE,self)
    return ret

  def arctan(self):
    ret = SXMatrix()
    ret.unary(ATAN_NODE,self)
    return ret

  def step(self):
    ret = SXMatrix()
    ret.unary(STEP_NODE,self)
    return ret

  def floor(self):
    ret = SXMatrix()
    ret.unary(FLOOR_NODE,self)
    return ret

  def ceil(self):
    ret = SXMatrix()
    ret.unary(CEIL_NODE,self)
    return ret

  def erf(self):
    ret = SXMatrix()
    ret.unary(ERF_NODE,self)
    return ret

  def fmin(self,other):
    ret = SXMatrix()
    ret.binary(FMIN_NODE,self,other)
    return ret

  def fmax(self,other):
    ret = SXMatrix()
    ret.binary(FMAX_NODE,self,other)
    return ret

  def __dealloc__(self):
    if self._sx_matrix_ref is not NULL:
      casadi_python.casadi_sx_matrix_delete(self._sx_matrix_ref)


cdef class FX:
  cpdef casadi_python.fx_ref _fx_ref

  def __cinit__(self, *args):
    self._fx_ref= casadi_python.casadi_fx_new()
    if self._fx_ref is NULL:
      python_exc.PyErr_NoMemory()
      
  def __dealloc__(self):
    if self._fx_ref is not NULL:
      casadi_python.casadi_fx_delete(self._fx_ref)

  cpdef print_string(self, STLString s):
    casadi_python.casadi_fx_print_string(self._fx_ref,s._string_ptr)

  def __repr__(self):
    s = STLString() 
    self.print_string(s)
    return str(s)

  def __str__(self):
    s = STLString() 
    self.print_string(s)
    return str(s)

  def set_input(self, np.ndarray[np.float64_t, ndim=1] v, *args):
    cdef int ind = 0
    cdef int order = 0
    if len(args) >= 1:
      ind = args[0]
    if len(args) >= 2:
      order = args[1]
    
    cdef int sz = -1  
    casadi_fx_input_size(self._fx_ref, ind, &sz)
    print sz
    assert sz == len(v), 'size mismatch'
    casadi_fx_setinput(self._fx_ref, ind, order, <double*>v.data)
    
  def set_output(self, np.ndarray[np.float64_t, ndim=1] v, *args):
    cdef int ind = 0
    cdef int order = 0
    if len(args) >= 1:
      ind = args[0]
    if len(args) >= 2:
      order = args[1]
    
    cdef int sz = -1  
    casadi_fx_output_size(self._fx_ref, ind, &sz)
    print sz
    assert sz == len(v), 'size mismatch'
    casadi_fx_setoutput(self._fx_ref, ind, order, <double*>v.data)
    
  def get_input(self, np.ndarray[np.float64_t, ndim=1] v, *args):
    cdef int ind = 0
    cdef int order = 0
    if len(args) >= 1:
      ind = args[0]
    if len(args) >= 2:
      order = args[1]
    
    cdef int sz = -1  
    casadi_fx_input_size(self._fx_ref, ind, &sz)
    assert sz == len(v), 'size mismatch'
    casadi_fx_getinput(self._fx_ref, ind, order, <double*>v.data)
    
  def get_output(self, np.ndarray[np.float64_t, ndim=1] v, *args):
    cdef int ind = 0
    cdef int order = 0
    if len(args) >= 1:
      ind = args[0]
    if len(args) >= 2:
      order = args[1]
    
    cdef int sz = -1  
    casadi_fx_output_size(self._fx_ref, ind, &sz)
    print sz
    assert sz == len(v), 'size mismatch'
    casadi_fx_getoutput(self._fx_ref, ind, order, <double*>v.data)

  def __setattr__num(self,name,np.ndarray[np.float64_t, ndim=1] value):
    casadi_fx_setoption_double(self._fx_ref,name,<double*>value.data,len(value))

  def __setattr__(self, name, value):
    assert isinstance(name,str)
    cdef double val

    if isinstance(value,str):
      casadi_fx_setoption_string(self._fx_ref,name,value)
    elif isinstance(value,int):
      val = value
      casadi_fx_setoption_double(self._fx_ref,name,&val,1)
    else:
      self.__setattr__num(name,value)

  cpdef getoption_string(self, str name, STLString s):
    casadi_fx_getoption_string(self._fx_ref,name,s._string_ptr)

  def __getattr__(self,name):
    s = STLString()
    self.getoption_string(name,s)
    return str(s)

  def print_options(self):
    casadi_fx_print_options(self._fx_ref)
    
  # initialize
  def init(self):
    casadi_fx_init(self._fx_ref)
     
   
    
    
    
    
cdef class MXVector:
  cpdef casadi_python.mx_vec _mx_vec

  def __cinit__(self, *args):
    self._mx_vec = casadi_python.casadi_mx_vec_new()
    if self._mx_vec is NULL:
      python_exc.PyErr_NoMemory()
      
  def __dealloc__(self):
    if self._mx_vec is not NULL:
      casadi_python.casadi_mx_vec_delete(self._mx_vec)

  cpdef push_back(self, MX ref):
    casadi_mx_vec_push_back(self._mx_vec, ref._mx_ref)
    
cdef class SXMatrixVector:
  cpdef casadi_python.sx_matrix_vec _sx_matrix_vec

  def __cinit__(self, *args):
    self._sx_matrix_vec = casadi_python.casadi_sx_matrix_vec_new()
    if self._sx_matrix_vec is NULL:
      python_exc.PyErr_NoMemory()
      
    # quick return if default constructor  
    if len(args)==0: return
    assert len(args)==1
    arg = args[0]

    print arg
    if isinstance(arg,list) or isinstance(arg,tuple):
      for x in arg:
        print "x = ", x
        print "SXMatrix(x) = ", SXMatrix(x)
        if isinstance(x,SXMatrix):
          self.push_back(x)
        else:
          self.push_back(SXMatrix(x))
    else:
      if isinstance(x,SXMatrix):
        self.push_back(x)
      else:
        self.push_back(SXMatrix(x))
        
      print "self = ", self
    
  def __dealloc__(self):
    if self._sx_matrix_vec is not NULL:
      casadi_python.casadi_sx_matrix_vec_delete(self._sx_matrix_vec)

  cpdef push_back(self, SXMatrix ref):
    casadi_sx_matrix_vec_push_back(self._sx_matrix_vec, ref._sx_matrix_ref)
    
  cpdef print_string(self, STLString s):
    casadi_python.casadi_sx_matrix_vec_print(self._sx_matrix_vec,s._string_ptr)

  def __repr__(self):
    s = STLString()
    self.print_string(s)
    return str(s)

  def __str__(self):
    s = STLString() 
    self.print_string(s)
    return str(s)

    

cdef class MXFunction(FX):
  cpdef mx_function(self, MXVector ip, MXVector op):
      casadi_python.casadi_mx_function(self._fx_ref,ip._mx_vec,op._mx_vec)
      
  def __cinit__(self, arg,res):
      #inputs
      ip = MXVector()
      if isinstance(arg,tuple):
        for el in arg:
          ip.push_back(el)
      else:
        ip.push_back(arg)

      #outputs
      op = MXVector()
      if isinstance(res,tuple):
        for el in res:
          op.push_back(el)
      else:
        op.push_back(res)

      self.mx_function(ip,op)


cdef class SXFunction(FX):
  cpdef sx_function(self, SXMatrixVector ip, SXMatrixVector op):
      casadi_python.casadi_sx_function(self._fx_ref,ip._sx_matrix_vec,op._sx_matrix_vec)
      
  def __cinit__(self, arg, res):
      if not isinstance(arg,SXMatrixVector): arg = SXMatrixVector(arg)
      if not isinstance(res,SXMatrixVector): res = SXMatrixVector(res)

      print "ok"
      print "arg = ", arg
      print "res = ", res

      ##inputs
      #ip = SXMatrixVector()
      #if isinstance(arg,tuple):
        #for el in arg:
          #ip.push_back(el)
      #else:
        #ip.push_back(arg)

      #outputs
      #op = SXMatrixVector()
      #if isinstance(res,tuple):
        #for el in res:
          #op.push_back(el)
      #else:
        #op.push_back(res)

      self.sx_function(arg,res)

cdef class IpoptSolver(FX):
  cpdef ipopt_solver(self, FX ffcn, FX gfcn, FX hfcn, FX jfcn):
      casadi_python.casadi_ipopt_solver(self._fx_ref,ffcn._fx_ref,gfcn._fx_ref,hfcn._fx_ref,jfcn._fx_ref)

  def __cinit__(self, ffcn, gfcn, hfcn, jfcn):
      self.ipopt_solver(ffcn,gfcn,hfcn,jfcn)

cdef class CVodesIntegrator(FX):
  cpdef cvodes_integrator(self, FX ffcn):
      casadi_cvodes_integrator(self._fx_ref,ffcn._fx_ref)

  def __cinit__(self, ffcn):
      self.cvodes_integrator(ffcn)

  def integrate(self, t_out):
      cdef double t = t_out
      casadi_integrator_integrate(self._fx_ref,t)

  def reset(self, with_sens):
      cdef int ws = with_sens
      casadi_integrator_reset(self._fx_ref,ws)
  
cdef vertcat(SXMatrix m, SXMatrixVector v):
  casadi_sx_matrix_vertcat(m._sx_matrix_ref, v._sx_matrix_vec)

cdef horzcat(SXMatrix m, SXMatrixVector v):
  casadi_sx_matrix_horzcat(m._sx_matrix_ref, v._sx_matrix_vec)

def concatenate(ars,axis=0):
  assert axis==0 or axis==1

  # Add arguments to an STL array
  v = SXMatrixVector()
  for x in ars:
    v.push_back(x)
    
  #Return matrix
  m = SXMatrix()
  if axis==0: vertcat(m,v)
  else:       horzcat(m,v)
  return m
