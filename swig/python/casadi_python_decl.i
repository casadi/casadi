%pythoncode %{

import contextlib

# copy / deepcopy support is added to the casadi value types via
# `%extend PrintableCommon` (their common base), not by shadowing the global
# `object` -- see further down in this file.

_swig_repr_default = _swig_repr
def _swig_repr(self):
  if hasattr(self,'repr'):
    return self.repr()
  else:
    return _swig_repr_default(self)

def DM_from_array(m, check_only=True):
  import numpy as np
  if isinstance(m, np.ndarray):
    if len(m.shape)>2:
      return False
    try:
      m = m.astype(float,casting="same_kind",copy=False)
    except:
      return False
    if check_only:
      return True
    else:
      shape = m.shape + (1, 1)
      nrow, ncol = shape[0], shape[1]
      return (nrow,ncol,m.flat)
  return False

def IM_from_array(m, check_only=True):
  import numpy as np
  if isinstance(m, np.ndarray):
    if len(m.shape)>2:
      return False
    try:
      m = m.astype(int,casting="same_kind",copy=False)
    except:
      return False
    if check_only:
      return True
    else:
      shape = m.shape + (1, 1)
      nrow, ncol = shape[0], shape[1]
      return (nrow,ncol,m.flat)
  return False

def SX_from_array(m, check_only=True):
  import numpy as np
  if isinstance(m, np.ndarray):
    if len(m.shape)>2:
      return False
    if m.dtype!=object: return None
    shape = m.shape + (1, 1)
    nrow, ncol = shape[0], shape[1]
    return (nrow,ncol,m.flat)
  return False

def DM_from_csc(m, check_only=True):
  if hasattr(m,"tocsc"):
    m = m.tocsc()
  if m.__class__.__name__ in ["csc_matrix", "csc_array"]:
    if len(m.shape)!=2: return False
    if check_only: return True
    return m.shape + (m.indptr.flat,m.indices.flat,m.data.flat)
  return False

%}

%feature("director:except") {
	if ($error != NULL) casadi::handle_director_exception();
}


%{
#define SWIG_FILE_WITH_INIT
%}


%typemap(in, doc="buffer(ro)", pystub_in="memoryview", noblock=1, fragment="casadi_all") (const double * a, casadi_int size) (Py_buffer _global_pybuf_ro) {
  if (PyObject_GetBuffer($input, &_global_pybuf_ro, PyBUF_SIMPLE) != 0) {
    SWIG_exception_fail(SWIG_TypeError, "Must supply a buffer-supporting object.");
  }
  $1 = static_cast<double*>(_global_pybuf_ro.buf);
  $2 = _global_pybuf_ro.len;
 }
%typemap(freearg) (const double * a, casadi_int size) {
  PyBuffer_Release(&_global_pybuf_ro);
 }

%typemap(in, doc="buffer(rw)", pystub_in="memoryview", noblock=1, fragment="casadi_all") (double * a, casadi_int size)  (Py_buffer _global_pybuf_rw) {
  if (PyObject_GetBuffer($input, &_global_pybuf_rw, PyBUF_WRITABLE) != 0) {
    SWIG_exception_fail(SWIG_TypeError, "Must supply a writable buffer-supporting object.");
  }
  $1 = static_cast<double*>(_global_pybuf_rw.buf);
  $2 = _global_pybuf_rw.len;
 }
%typemap(freearg) (double * a, casadi_int size) {
  PyBuffer_Release(&_global_pybuf_rw);
 }

// Directorin typemap; as output
%typemap(directorin, noblock=1, fragment="casadi_all") (const double** arg, const std::vector<casadi_int>& sizes_arg) (PyObject* my_tuple) {
  PyObject * arg_tuple = PyTuple_New($2.size());
  for (casadi_int i=0;i<$2.size();++i) {

    PyObject* buf = $1[i] ? PyMemoryView_FromMemory(reinterpret_cast<char*>(const_cast<double*>($1[i])), $2[i]*sizeof(double), PyBUF_READ) : SWIG_Py_Void();
    PyTuple_SET_ITEM(arg_tuple, i, buf);
  }
  $input = arg_tuple;
}

%typemap(directorin, noblock=1, fragment="casadi_all") (double** res, const std::vector<casadi_int>& sizes_res) {
  PyObject* res_tuple = PyTuple_New($2.size());
  for (casadi_int i=0;i<$2.size();++i) {
    PyObject* buf = $1[i] ? PyMemoryView_FromMemory(reinterpret_cast<char*>(const_cast<double*>($1[i])), $2[i]*sizeof(double), PyBUF_WRITE) : SWIG_Py_Void();
    PyTuple_SET_ITEM(res_tuple, i, buf);
  }
  $input = res_tuple;
}

%typemap(in, doc="void*", pystub_in="Any", noblock=1, fragment="casadi_all") void* raw {
  $1 = PyCapsule_GetPointer($input, NULL);
}

%typemap(out, doc="void*", pystub_out="Any", noblock=1, fragment="casadi_all") void* {
  $result = PyCapsule_New($1, NULL,NULL);
}

%stubcode %{class ArrayInterface(Generic[_T]):
    """numpy-semantics array view over a casadi DM/SX/MX.  Constructing one
    dispatches to the backing-typed subclass ArrayInterfaceDM/SX/MX."""
    __array_priority__: float
    @overload
    def __new__(cls, value: SX, ndim: int | None = None) -> ArrayInterfaceSX: ...
    @overload
    def __new__(cls, value: MX, ndim: int | None = None) -> ArrayInterfaceMX: ...
    @overload
    def __new__(cls, value: object = None, ndim: int | None = None) -> ArrayInterfaceDM: ...
    def __init__(self, value: object, ndim: int | None = None) -> None: ...
    @property
    def shape(self) -> tuple[int, ...]: ...
    @property
    def ndim(self) -> int: ...
    @property
    def size(self) -> int: ...
    def to_casadi(self) -> _T: ...
    @property
    def T(self) -> ArrayInterface[_T]: ...
    def reshape(self, *shape: int) -> ArrayInterface[_T]: ...
    def transpose(self, *axes: int) -> ArrayInterface[_T]: ...
    def squeeze(self, axis: int | tuple[int, ...] | None = None) -> ArrayInterface[_T]: ...
    def flatten(self) -> ArrayInterface[_T]: ...
    def ravel(self) -> ArrayInterface[_T]: ...
    @property
    def flat(self) -> Any: ...
    def to_DM(self) -> DM: ...
    def sum(self, axis: int | tuple[int, ...] | None = None) -> ArrayInterface[_T]: ...
    def mean(self, axis: int | tuple[int, ...] | None = None) -> ArrayInterface[_T]: ...
    def __getitem__(self, idx: Any) -> ArrayInterface[_T]: ...
    def __setitem__(self, idx: Any, value: Any) -> None: ...
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[ArrayInterface[_T]]: ...
    def __array__(self, dtype: Any = None, copy: bool | None = None) -> NDArray[Any]: ...
    def __array_ufunc__(self, ufunc: Any, method: builtins.str, *inputs: Any, **kwargs: Any) -> Any: ...
    def __array_function__(self, func: Any, types: Any, args: Any, kwargs: Any) -> Any: ...
    def __float__(self) -> float: ...
    def __int__(self) -> int: ...
    def __bool__(self) -> bool: ...
    def __add__(self, o: Any) -> ArrayInterface[_T]: ...
    def __radd__(self, o: Any) -> ArrayInterface[_T]: ...
    def __sub__(self, o: Any) -> ArrayInterface[_T]: ...
    def __rsub__(self, o: Any) -> ArrayInterface[_T]: ...
    def __mul__(self, o: Any) -> ArrayInterface[_T]: ...
    def __rmul__(self, o: Any) -> ArrayInterface[_T]: ...
    def __truediv__(self, o: Any) -> ArrayInterface[_T]: ...
    def __rtruediv__(self, o: Any) -> ArrayInterface[_T]: ...
    def __floordiv__(self, o: Any) -> ArrayInterface[_T]: ...
    def __rfloordiv__(self, o: Any) -> ArrayInterface[_T]: ...
    def __mod__(self, o: Any) -> ArrayInterface[_T]: ...
    def __rmod__(self, o: Any) -> ArrayInterface[_T]: ...
    def __pow__(self, o: Any) -> ArrayInterface[_T]: ...
    def __rpow__(self, o: Any) -> ArrayInterface[_T]: ...
    def __matmul__(self, o: Any) -> ArrayInterface[_T]: ...
    def __rmatmul__(self, o: Any) -> ArrayInterface[_T]: ...
    def __lt__(self, o: Any) -> ArrayInterface[_T]: ...
    def __le__(self, o: Any) -> ArrayInterface[_T]: ...
    def __gt__(self, o: Any) -> ArrayInterface[_T]: ...
    def __ge__(self, o: Any) -> ArrayInterface[_T]: ...
    def __eq__(self, o: object) -> ArrayInterface[_T]: ...  # type: ignore[override]
    def __ne__(self, o: object) -> ArrayInterface[_T]: ...  # type: ignore[override]
    __hash__: None  # type: ignore[assignment]
    def __neg__(self) -> ArrayInterface[_T]: ...
    def __pos__(self) -> ArrayInterface[_T]: ...
    def __abs__(self) -> ArrayInterface[_T]: ...
    def prod(self, axis: int | tuple[int, ...] | None = None) -> ArrayInterface[_T]: ...
    def max(self, axis: int | tuple[int, ...] | None = None) -> ArrayInterface[_T]: ...
    def min(self, axis: int | tuple[int, ...] | None = None) -> ArrayInterface[_T]: ...
    def ptp(self, axis: int | tuple[int, ...] | None = None) -> ArrayInterface[_T]: ...
    def all(self, axis: int | tuple[int, ...] | None = None) -> ArrayInterface[_T]: ...
    def any(self, axis: int | tuple[int, ...] | None = None) -> ArrayInterface[_T]: ...
    def var(self, axis: int | tuple[int, ...] | None = None, ddof: int = 0) -> ArrayInterface[_T]: ...
    def std(self, axis: int | tuple[int, ...] | None = None, ddof: int = 0) -> ArrayInterface[_T]: ...
    @classmethod
    def DM(cls, value: Any, shape: Any = ...) -> ArrayInterfaceDM: ...
    @classmethod
    def SX(cls, value: Any, shape: Any = ...) -> ArrayInterfaceSX: ...
    @classmethod
    def MX(cls, value: Any, shape: Any = ...) -> ArrayInterfaceMX: ...

class ArrayInterfaceDM(ArrayInterface[DM]):
    def __DM__(self) -> DM: ...

class ArrayInterfaceSX(ArrayInterface[SX]):
    def __SX__(self) -> SX: ...

class ArrayInterfaceMX(ArrayInterface[MX]):
    def __MX__(self) -> MX: ...

%}

%pythoncode %{
if __name__ != "casadi.casadi":
  raise Exception("""
            CasADi is not running from its package context.

            You probably specified the wrong casadi directory.

            When setting PYTHONPATH or sys.path.append,
            take care not to add a trailing '/casadi'.

        """)

def swigtypeconvertor(*args):
  return swig_typename_convertor_python2cpp(args)

def swig_typename_convertor_python2cpp(a):
  try:
    import numpy as np
  except:
    class NoExist:
      pass
    class Temp(object):
      ndarray = NoExist
    np = Temp()
  if isinstance(a,list):
    if len(a)>0:
      return "[%s]" % "|".join(set([swig_typename_convertor_python2cpp(i) for i in a]))
    else:
      return "[]"
  elif isinstance(a,tuple):
    return "(%s)" % ",".join([swig_typename_convertor_python2cpp(i) for i in a])
  elif isinstance(a,np.ndarray):
    return "np.array(%s)" % ",".join(set([swig_typename_convertor_python2cpp(i) for i in np.array(a).flatten().tolist()]))
  elif isinstance(a,dict):
    if len(a)>0:
      return "|".join(set([swig_typename_convertor_python2cpp(i) for i in a.keys()])) +":"+ "|".join(set([swig_typename_convertor_python2cpp(i) for i in a.values()]))
    else:
      return "dict"
  return type(a).__name__
%}

#ifdef WITH_PYTHON_INTERRUPTS
%{
#include <pythonrun.h>
void SigIntHandler(casadi_int) {
  std::cerr << "Keyboard Interrupt" << std::endl;
  signal(SIGINT, SIG_DFL);
  kill(getpid(), SIGINT);
}
%}

%init %{
PyOS_setsig(SIGINT, SigIntHandler);
%}
#endif // WITH_PYTHON_INTERRUPTS

%pythoncode "numpy_bridge.py"

%pythoncode "casadi_nparray.py"

// No inf/pi stubs: they would clash with numpy's Final declarations

// Stubs for the numpy-flavoured aliases above
%define %stub_unary_math(NAME)
%stub_overload_func(NAME, DM, x: _DM)
%stub_overload_func(NAME, SX, x: _SX)
%stub_overload_func(NAME, MX, x: _MX)
%enddef
%stub_unary_math(arcsin)
%stub_unary_math(arccos)
%stub_unary_math(arctan)
%stub_unary_math(arctanh)
%stub_unary_math(arcsinh)
%stub_unary_math(arccosh)
/* Re-exported from numpy at runtime (`from numpy import pi, inf`). */
#ifdef SWIG_STUBS_ENABLED
%stubcode %{pi: float
inf: float
%}
%stub_overload_func(arctan2, DM, x: _DM, y: _DM)
%stub_overload_func(arctan2, SX, x: _SX, y: _SX)
%stub_overload_func(arctan2, MX, x: _MX, y: _MX)
#endif // SWIG_STUBS_ENABLED

%ignore T;

%rename(logic_and) casadi_and;
%rename(logic_or) casadi_or;
%rename(logic_not) casadi_not;
%rename(logic_all) casadi_all;
%rename(logic_any) casadi_any;
%rename(fabs) casadi_abs;

// Concatenations
%rename(_veccat) casadi_veccat;
%rename(_vertcat) casadi_vertcat;
%rename(_horzcat) casadi_horzcat;
%rename(_diagcat) casadi_diagcat;
%pythoncode %{
def veccat(*args):
    try:
        if len(args)==0:
            return DM(0,1)
    except:
        pass
    return _veccat(args)
def vertcat(*args):
    try:
        if len(args)==0:
            return DM(0,1)
    except:
        pass
    return _vertcat(args)
def horzcat(*args):
    try:
        if len(args)==0:
            return DM(1,0)
    except:
        pass
    return _horzcat(args)
def diagcat(*args):
    try:
        if len(args)==0:
            return DM(0,0)
    except:
        pass
    return _diagcat(args)
def vvcat(args):
    try:
        if len(args)==0:
            return DM(0,1)
    except:
        pass
    return _veccat(args)
def vcat(args):
    try:
        if len(args)==0:
            return DM(0,1)
    except:
        pass
    return _vertcat(args)
def hcat(args):
    try:
        if len(args)==0:
            return DM(1,0)
    except:
        pass
    return _horzcat(args)
def dcat(args):
    try:
        if len(args)==0:
            return DM(0,0)
    except:
        pass
    return _diagcat(args)
%}

// Stubs for the *cat wrappers above, narrowest first
%stub_overload_func(veccat,  Sparsity, *args: Sparsity)
%stub_overload_func(veccat,  DM,       *args: _DM)
%stub_overload_func(veccat,  SX,       *args: _SX)
%stub_overload_func(veccat,  MX,       *args: _MX)
%stub_overload_func(vertcat, Sparsity, *args: Sparsity)
%stub_overload_func(vertcat, DM,       *args: _DM)
%stub_overload_func(vertcat, SX,       *args: _SX)
%stub_overload_func(vertcat, MX,       *args: _MX)
%stub_overload_func(horzcat, Sparsity, *args: Sparsity)
%stub_overload_func(horzcat, DM,       *args: _DM)
%stub_overload_func(horzcat, SX,       *args: _SX)
%stub_overload_func(horzcat, MX,       *args: _MX)
%stub_overload_func(diagcat, Sparsity, *args: Sparsity)
%stub_overload_func(diagcat, DM,       *args: _DM)
%stub_overload_func(diagcat, SX,       *args: _SX)
%stub_overload_func(diagcat, MX,       *args: _MX)
%stub_overload_func(vvcat, Sparsity, args: Sequence[Sparsity])
%stub_overload_func(vvcat, DM,       args: Sequence[_DM])
%stub_overload_func(vvcat, SX,       args: Sequence[_SX])
%stub_overload_func(vvcat, MX,       args: Sequence[_MX])
%stub_overload_func(vcat,  Sparsity, args: Sequence[Sparsity])
%stub_overload_func(vcat,  DM,       args: Sequence[_DM])
%stub_overload_func(vcat,  SX,       args: Sequence[_SX])
%stub_overload_func(vcat,  MX,       args: Sequence[_MX])
%stub_overload_func(hcat,  Sparsity, args: Sequence[Sparsity])
%stub_overload_func(hcat,  DM,       args: Sequence[_DM])
%stub_overload_func(hcat,  SX,       args: Sequence[_SX])
%stub_overload_func(hcat,  MX,       args: Sequence[_MX])
%stub_overload_func(dcat,  Sparsity, args: Sequence[Sparsity])
%stub_overload_func(dcat,  DM,       args: Sequence[_DM])
%stub_overload_func(dcat,  SX,       args: Sequence[_SX])
%stub_overload_func(dcat,  MX,       args: Sequence[_MX])

// Non-fatal errors (returning NotImplemented singleton)
%feature("python:maybecall") casadi_plus;
%feature("python:maybecall") casadi_minus;
%feature("python:maybecall") casadi_times;
%feature("python:maybecall") casadi_rdivide;
%feature("python:maybecall") casadi_lt;
%feature("python:maybecall") casadi_le;
%feature("python:maybecall") casadi_eq;
%feature("python:maybecall") casadi_ne;
%feature("python:maybecall") casadi_power;
%feature("python:maybecall") casadi_atan2;
%feature("python:maybecall") casadi_and;
%feature("python:maybecall") casadi_or;
%feature("python:maybecall") casadi_mod;
%feature("python:maybecall") casadi_copysign;
%feature("python:maybecall") casadi_constpow;

%ignore casadi_mod;
%rename(__bool__) __nonzero__;

%pythoncode %{
class NZproxy:
  def __init__(self,matrix):
    self.matrix = matrix

  def __getitem__(self,s):
    return self.matrix.get_nz(False, s)

  def __setitem__(self,s,val):
    return self.matrix.set_nz(val, False, s)

  def __len__(self):
    return self.matrix.nnz()

  def __iter__(self):
    for i in range(len(self)):
      yield self[i]

%}

// NZproxy (%pythoncode), generic in the owning matrix type
#ifdef SWIG_STUBS_ENABLED
%stubcode %{class NZproxy(Generic[_TNZ]):
%}
%stub_method(__init__,    None, matrix: _TNZ)
%stub_overload_method_selftyped(__getitem__, MX, NZproxy[MX], s: _MIndex | MX)
%stub_overload_method(__getitem__, _TNZ,  s: _MIndex)
%stub_overload_method_selftyped(__setitem__, None, NZproxy[MX], s: _MIndex | MX, val: float | MX | Sequence[float])
%stub_overload_method(__setitem__, None, s: _MIndex, val: float | _TNZ | Sequence[float])
%stub_method0(__len__,  int)
%stub_method0(__iter__, Iterator[_T])
#endif

%define %python_array_wrappers(arraypriority, arraydtype)
%pythoncode %{

  __array_priority__ = arraypriority

  def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
      return _numpy_ufunc_dispatch(self, ufunc, method, inputs, kwargs)

  def __array_function__(self, func, types, args, kwargs):
      return _numpy_array_function_dispatch(self, func, types, args, kwargs)

  def __array__(self, *args, **kwargs):
      import numpy as _n
      try:
          arr = self.full()
      except Exception:
          if self.is_scalar(True):
              # Box symbolic scalars in an object array (#2743).
              E = _n.empty((), dtype=object)
              E[()] = self
              return E
          raise TypeError(
              "Implicit conversion of symbolic CasADi type to numeric matrix not supported.\n"
              "This may occur when you pass a CasADi object to a numpy function.\n"
              "Use an equivalent CasADi function instead of that numpy function.")
      dtype = kwargs.get("dtype")
      if dtype is not None and dtype is not _n.double:
          return _n.array(arr, dtype=dtype)
      return arr

%}
// __array__ makes DM/SX/MX numpy ArrayLike
%stub_method(__array__, %arg(NDArray[arraydtype]), *args: Any, **kwargs: Any)
#ifdef SWIG_STUBS_ENABLED
%stubcode %{    __array_priority__: float
%}
#endif // SWIG_STUBS_ENABLED
%stub_method(__array_ufunc__, Any, ufunc: Any, method: builtins.str, *inputs: Any, **kwargs: Any)
%stub_method(__array_function__, Any, func: Any, types: Any, args: Any, kwargs: Any)
%enddef

%feature("python:stub_skip") casadi::Callback::eval;
%feature("python:stub_skip") casadi::Callback::eval_buffer;
%feature("python:stub_skip") casadi::Callback::get_sparsity_in;
%feature("python:stub_skip") casadi::Callback::get_sparsity_out;

%rename("_pop_type") casadi::DeserializerBase::pop_type;
%rename("%(regex:/(SERIALIZED_\w+)/_\\1/)s", regextarget=1, fullname=1) "casadi::SerializerBase::SERIALIZED_\w+";

  %rename(_variable) casadi::Opti::variable;
  %rename(_parameter) casadi::Opti::parameter;
  %rename(_subject_to) casadi::Opti::subject_to;
  %extend casadi::Opti {
    %pythoncode %{
      def parameter(self,*args):
        import sys
        import os
        try:
            frame = sys._getframe(1)
        except:
            frame = {}
        meta = {} if frame is None else {"stacktrace": [{"file":os.path.abspath(frame.f_code.co_filename),"line":frame.f_lineno,"name":frame.f_code.co_name}]}
        ret = self._parameter(*args)
        if len(meta)>0:
            self.update_user_dict(ret, meta)
        return ret

      def variable(self,*args):
        import sys
        import os
        try:
            frame = sys._getframe(1)
        except:
            frame = {}
        meta = {} if frame is None else {"stacktrace": [{"file":os.path.abspath(frame.f_code.co_filename),"line":frame.f_lineno,"name":frame.f_code.co_name}]}
        ret = self._variable(*args)
        if len(meta)>0:
            self.update_user_dict(ret, meta)
        return ret

      def subject_to(self,*args):
        if len(args)==0:
          return self._subject_to()
        import sys
        import os
        stacktrace = []
        for i in range(1,10000):
          try:
            frame = sys._getframe(i)
            stacktrace.append({"file":os.path.abspath(frame.f_code.co_filename),"line":frame.f_lineno,"name":frame.f_code.co_name})
          except Exception as e:
            break
        args = list(args)
        if len(args)==3 and isinstance(args[2],dict):
          args[2] = dict(args[2])
          if "stacktrace" not in args[2]:
            args[2]["stacktrace"] = stacktrace
        elif len(args)==2 and isinstance(args[1],dict):
          args[1] = dict(args[1])
          if "stacktrace" not in args[1]:
            args[1]["stacktrace"] = stacktrace
        elif len(args)==1:
          args = [args[0], {"stacktrace": stacktrace}]
        elif len(args)==2:
          args = [args[0], args[1], {"stacktrace": stacktrace}]
        ret = self._subject_to(*args)
        return ret
    %}
    /* Stubs for the three %pythoncode methods above. */
    %stub_overload_method(variable, MX, n: int = ..., m: int = ..., attribute: builtins.str = ...)
    %stub_overload_method(variable, MX, sp: Sparsity | DM, attribute: builtins.str = ...)
    %stub_overload_method(variable, MX, symbol: MX | SX | DM, attribute: builtins.str = ...)
    %stub_overload_method(parameter, MX, n: int = ..., m: int = ..., attribute: builtins.str = ...)
    %stub_overload_method(parameter, MX, sp: Sparsity | DM, attribute: builtins.str = ...)
    %stub_overload_method(parameter, MX, symbol: MX | SX | DM, attribute: builtins.str = ...)
    // subject_to(g[, linear_scale][, options])
    %stub_overload_method0(subject_to, None)
    %stub_overload_method(subject_to, None, g: MX | SX | DM | float)
    %stub_overload_method(subject_to, None, g: Sequence[MX | SX | DM | float])
    %stub_overload_method(subject_to, None, g: MX | SX | DM | float, options: Mapping[builtins.str, _GenericType])
    %stub_overload_method(subject_to, None, g: Sequence[MX | SX | DM | float], options: Mapping[builtins.str, _GenericType])
    %stub_overload_method(subject_to, None, g: MX | SX | DM | float, linear_scale: _DM)
    %stub_overload_method(subject_to, None, g: Sequence[MX | SX | DM | float], linear_scale: _DM)
    %stub_overload_method(subject_to, None, g: MX | SX | DM | float, linear_scale: _DM, options: Mapping[builtins.str, _GenericType])
    %stub_overload_method(subject_to, None, g: Sequence[MX | SX | DM | float], linear_scale: _DM, options: Mapping[builtins.str, _GenericType])
  }

%extend casadi::Sparsity{
    %pythoncode %{
        @property
        def shape(self):
            return (self.size1(),self.size2())

        @property
        def T(self):
            return _casadi.transpose(self)

        def __array__(self,*args,**kwargs):
            return DM.ones(self).full()
    %}
};

%pythoncode %{

try:
  import numpy

  def constpow(x,y):
    pass

  constpow=numpy.frompyfunc(constpow,2,1)
except:
  pass
%}

make_property(casadi::Opti, casadi_solver, Function);

// Explicit conversion to double and casadi_int
%rename(__float__) operator double;
%rename(__int__) operator casadi_int;

// full() and sparse() call the Python API: keep the GIL
%feature("nothread") casadi::Matrix<double>::full;
%feature("nothread") casadi::Matrix<double>::sparse;
