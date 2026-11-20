namespace casadi{
%extend Sparsity {
  %pythoncode %{
    def __setstate__(self, state):
        self.__init__(Sparsity.deserialize(state["serialization"]))

    def __getstate__(self):
        return {"serialization": self.serialize()}
  %}
  %stub_method0(__getstate__, %arg(dict[builtins.str, Any]))
  %stub_method(__setstate__, None, state: Mapping[builtins.str, Any])
  %stub_method(__array__, %arg(NDArray[Any]), *args: Any, **kwargs: Any)
  // Getter-backed attributes without a typemap doc
  %stub_attr(T, Sparsity)
  %stub_attr(shape, tuple[int, int]) // TYPE is variadic; inner comma is fine here
}
%extend Matrix<SXElem> {
  %pythoncode %{
    def __setstate__(self, state):
      ctx = _current_unpickle_context()
      if not ctx:
        raise Exception("Cannot unpickle SX objects without a casadi context. " +
          "Use something like:\n"+
          "with ca.global_unpickle_context(): \n"+
          "  f_ref = pickle.load(open(filename,'rb'))")
      ctx.decode(state)
      self.__init__(ctx.unpack())

    def __getstate__(self):
      ctx = _current_pickle_context()
      if not ctx:
        raise Exception("Cannot pickle SX objects without a casadi context. " +
          "Use something like:\n"+
          "with ca.global_pickle_context(): \n"+
          "  pickle.dump(f,open(filename,'wb'))")
      ctx.pack(self)
      return ctx.encode()
  %}
  %stub_method0(__getstate__, Any)
  %stub_method(__setstate__, None, state: Any)
}
%extend MX {
  %pythoncode %{
    def __setstate__(self, state):
      ctx = _current_unpickle_context()
      if not ctx:
        raise Exception("Cannot unpickle MX objects without a casadi context. " +
          "Use something like:\n"+
          "with ca.global_unpickle_context(): \n"+
          "  f_ref = pickle.load(open(filename,'rb'))")
      ctx.decode(state)
      self.__init__(ctx.unpack())

    def __getstate__(self):
      ctx = _current_pickle_context()
      if not ctx:
        raise Exception("Cannot pickle MX objects without a casadi context. " +
          "Use something like:\n"+
          "with ca.global_pickle_context(): \n"+
          "  pickle.dump(f,open(filename,'wb'))")
      ctx.pack(self)
      return ctx.encode()
  %}
  %stub_method0(__getstate__, Any)
  %stub_method(__setstate__, None, state: Any)
}

} // namespace casadi

namespace casadi{
%extend Matrix<double> {

%python_array_wrappers(999.0, np.float64)

%pythoncode %{
  def tocsc(self):
    import numpy as np
    import warnings
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      from scipy.sparse import csc_matrix
    return csc_matrix( (self.nonzeros(),self.row(),self.colind()), shape = self.shape, dtype=np.double )
  def tocsc_array(self):
    import numpy as np
    import warnings
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      from scipy.sparse import csc_array
    return csc_array( (self.nonzeros(),self.row(),self.colind()), shape = self.shape, dtype=np.double )
  def toarray(self,simplify=False):
    import numpy as np
    if simplify:
      if self.is_scalar():
        return float(self)
      elif self.is_vector():
        return np.array(self.elements())
    # elements() is column-major; avoids materializing the transpose
    return np.ascontiguousarray(np.array(self.elements()).reshape(self.shape, order='F'))
  def to_numpy(self):
    # Pandas-style densification hook.  matplotlib's
    # cbook._unpack_to_numpy looks for `to_numpy()` before falling back
    # to iteration, so defining it here makes plt.plot(DM) work without
    # forcing densification through NEP-18 shape ops like np.atleast_1d.
    return self.full()
%}
// scipy.sparse result, left untyped
%stub_method0(tocsc,    Any)
%stub_method0(tocsc_array, Any)
%stub_method(toarray,   %arg(NDArray[np.float64] | float), simplify: bool = ...)
%stub_method0(to_numpy, %arg(NDArray[np.float64]))


%pythoncode %{
  def __bool__(self):
    if self.numel()!=1:
      raise Exception("Only a scalar can be cast to a float")
    if self.nnz()==0:
      return False
    return float(self)!=0
%}

}; // extend Matrix<double>


// Logic for pickling

%extend Matrix<double> {

  %pythoncode %{
    def __setstate__(self, state):
        self.__init__(DM.deserialize(state["serialization"]))

    def __getstate__(self):
        return {"serialization": self.serialize()}
  %}
  %stub_method0(__getstate__, %arg(dict[builtins.str, Any]))
  %stub_method(__setstate__, None, state: Mapping[builtins.str, Any])

}


%extend Function {

  %pythoncode %{
    def __setstate__(self, state):
        self.__init__(Function.deserialize(state["serialization"]))

    def __getstate__(self):
        return {"serialization": self.serialize()}
  %}
  %stub_method0(__getstate__, %arg(dict[builtins.str, Any]))
  %stub_method(__setstate__, None, state: Mapping[builtins.str, Any])

}


} // namespace casadi

#include <arrayobject.h>
%template()    std::vector<PyObject*>;

%pythoncode %{
def attach_return_type(f,t):
  if not(hasattr(f,'func_annotations')):
    f.func_annotations = {}
  if not(isinstance(getattr(f,'func_annotations'),dict)):
    raise Exception("Cannot annotate this python Method to be a sparsitygenerator. Method has func_annotations attribute with unknown type.")
  f.func_annotations["return"] = t
  return f

def pyevaluate(f):
  return attach_return_type(f,None)

def pycallback(f):
  return attach_return_type(f,int)


def pyfunction(inputs,outputs):
  def wrap(f):

    @pyevaluate
    def fcustom(f2):
      res = f([f2.getInput(i) for i in range(f2.n_in())])
      if not isinstance(res,list):
        res = [res]
      for i in range(f2.n_out()):
        f2.setOutput(res[i],i)
    import warnings

    with warnings.catch_warnings():
      warnings.filterwarnings("ignore",category=DeprecationWarning)
      Fun = CustomFunction("CustomFunction",fcustom,inputs,outputs)
      return Fun

  return wrap

def PyFunction(name, obj, inputs, outputs, opts={}):
    @pyevaluate
    def fcustom(f):
      res = [f.getOutput(i) for i in range(f.n_out())]
      obj.evaluate([f.getInput(i) for i in range(f.n_in())],res)
      for i in range(f.n_out()): f.setOutput(res[i], i)

    import warnings

    with warnings.catch_warnings():
      warnings.filterwarnings("ignore",category=DeprecationWarning)
      return CustomFunction("CustomFunction", fcustom,
                            inputs, outputs, opts)

%}

namespace casadi{
%extend Function {
  %pythoncode %{
    def __call__(self, *args, **kwargs):
      # Either named inputs or ordered inputs
      if len(args)>0 and len(kwargs)>0:
        raise SyntaxError('Function evaluation requires all arguments to be named or none')
      if len(args)>0:
        # Ordered inputs -> return tuple
        ret = self.call(args)
        if len(ret)==0:
          return None
        elif len(ret)==1:
          return ret[0]
        else:
          return tuple(ret)
      else:
        # Named inputs -> return dictionary
        return self.call(kwargs)

    def buffer(self):
      """
      Create a FunctionBuffer object for evaluating with minimal overhead

      """
      import functools
      fb = FunctionBuffer(self)
      caller = functools.partial(_casadi._function_buffer_eval, fb._self())
      return (fb, caller)
  %}
  // buffer() returns (FunctionBuffer, callable)
  %stub_method0(buffer, %arg(tuple[FunctionBuffer, Any]))
  // __call__ narrows per homogeneous input type; f() and f(**kw) give dicts
  %stub_overload_method(__call__, DM, __arg0: _DM, *args: _DM)
  %stub_overload_method(__call__, SX, __arg0: _SX, *args: _SX)
  %stub_overload_method(__call__, MX, __arg0: _MX, *args: _MX)
  %stub_overload_method(__call__, %arg(dict[builtins.str, DM]), **kwargs: _DM)
  %stub_overload_method(__call__, %arg(dict[builtins.str, SX]), **kwargs: _SX)
  %stub_overload_method(__call__, %arg(dict[builtins.str, MX]), **kwargs: _MX)
  %stub_overload_method(__call__, %arg(dict[builtins.str, DM | SX | MX]), **kwargs: float | DM | SX | MX)
 }

}

%pythoncode %{
def _array_priority_yield(x, y):
    # numpy-style operator deferral: a binary op on a casadi value yields
    # (returns NotImplemented) to a NON-casadi operand that advertises a
    # higher __array_priority__ -- e.g. casadi.ArrayInterface -- so its
    # reflected operator runs and the higher-priority type controls the
    # result.  casadi's own DM/SX/MX promotion is untouched: those are
    # excluded, and their mixing is handled natively by the _casadi.* call.
    py = getattr(y, "__array_priority__", None)
    if py is None or py <= getattr(x, "__array_priority__", 0.0):
        return False
    return not isinstance(y, (DM, SX, MX, Sparsity))
%}
namespace casadi {
  %extend GenericExpressionCommon {
// Upstream SWIG raises TypeError on a failed overload dispatch, the fork NotImplementedError
#ifdef SWIG_STUBS_ENABLED
    %pythoncode %{
      def __hash__(self):
        try:
          return self.element_hash()
        except:
          return SharedObject.__hash__(self)
      def __matmul__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        try:
          return _casadi.mtimes(x, y)
        except NotImplementedError:
          return NotImplemented
      def __rmatmul__(x, y):
        try:
          return _casadi.mtimes(y, x)
        except NotImplementedError:
          return NotImplemented
      def __imatmul__(x, y):
        try:
          return _casadi.mtimes(x, y)
        except NotImplementedError:
          return NotImplemented
    %}
#else
    %pythoncode %{
      def __hash__(self):
        try:
          return self.element_hash()
        except:
          return SharedObject.__hash__(self)
      def __matmul__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        try:
          return _casadi.mtimes(x, y)
        except (NotImplementedError, TypeError):
          return NotImplemented
      def __rmatmul__(x, y):
        try:
          return _casadi.mtimes(y, x)
        except (NotImplementedError, TypeError):
          return NotImplemented
      def __imatmul__(x, y):
        try:
          return _casadi.mtimes(x, y)
        except (NotImplementedError, TypeError):
          return NotImplemented
    %}
#endif
    %stub_method0(__hash__, int)
    %stub_CasadiMatrix_binop(__matmul__)
    %stub_CasadiMatrix_binop(__rmatmul__)
    %stub_CasadiMatrix_binop(__imatmul__)
  }
}
namespace casadi {
  %extend GenericExpressionCommon {
    %pythoncode %{
      def __add__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        return _casadi.plus(x, y)
      def __radd__(x, y): return _casadi.plus(y, x)
      def __sub__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        return _casadi.minus(x, y)
      def __rsub__(x, y): return _casadi.minus(y, x)
      def __mul__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        return _casadi.times(x, y)
      def __rmul__(x, y): return _casadi.times(y, x)
      def __truediv__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        return _casadi.rdivide(x, y)
      def __rtruediv__(x, y): return _casadi.rdivide(y, x)
      def __floordiv__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        return _casadi.floor(_casadi.rdivide(x, y))
      def __rfloordiv__(x, y): return _casadi.floor(_casadi.rdivide(y, x))
      def __mod__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        return x - y * _casadi.floor(_casadi.rdivide(x, y))
      def __rmod__(x, y):
        return y - x * _casadi.floor(_casadi.rdivide(y, x))
      def __divmod__(x, y):
        q = _casadi.floor(_casadi.rdivide(x, y))
        return (q, x - y * q)
      def __rdivmod__(x, y):
        q = _casadi.floor(_casadi.rdivide(y, x))
        return (q, y - x * q)
      def __lt__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        return _casadi.lt(x, y)
      def __rlt__(x, y): return _casadi.lt(y, x)
      def __le__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        return _casadi.le(x, y)
      def __rle__(x, y): return _casadi.le(y, x)
      def __gt__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        return _casadi.lt(y, x)
      def __rgt__(x, y): return _casadi.lt(x, y)
      def __ge__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        return _casadi.le(y, x)
      def __rge__(x, y): return _casadi.le(x, y)
      def __eq__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        r = _casadi.eq(x, y)
        if r is NotImplemented and isinstance(x, SX) and isinstance(y, MX):
          raise Exception("Cannot compare SX and MX objects for equality")
        return r
      def __ne__(x, y):
        if _array_priority_yield(x, y): return NotImplemented
        r = _casadi.ne(x, y)
        if r is NotImplemented and isinstance(x, SX) and isinstance(y, MX):
          raise Exception("Cannot compare SX and MX objects for inequality")
        return r
      def __req__(x, y):
        r = _casadi.eq(y, x)
        if r is NotImplemented and isinstance(x, SX) and isinstance(y, MX):
          raise Exception("Cannot compare SX and MX objects for equality")
        return r
      def __rne__(x, y):
        r = _casadi.ne(y, x)
        if r is NotImplemented and isinstance(x, SX) and isinstance(y, MX):
          raise Exception("Cannot compare SX and MX objects for inequality")
        return r
      def __pow__(x, n, modulo=None):
        if modulo is None and _array_priority_yield(x, n): return NotImplemented
        p = _casadi.power(x, n)
        if modulo is None:
          return p
        return p - modulo * _casadi.floor(_casadi.rdivide(p, modulo))
      def __rpow__(n, x): return _casadi.power(x, n)
      def __round__(x, ndigits=None):
        if ndigits is None:
          return _casadi.if_else(x >= 0, _casadi.floor(x + 0.5),
                                 _casadi.ceil(x - 0.5))
        s = 10.0 ** ndigits
        y = x * s
        return _casadi.if_else(y >= 0, _casadi.floor(y + 0.5),
                               _casadi.ceil(y - 0.5)) / s
      def __trunc__(x):
        return _casadi.if_else(x >= 0, _casadi.floor(x), _casadi.ceil(x))
      def __floor__(x): return _casadi.floor(x)
      def __ceil__(x):  return _casadi.ceil(x)
      def __abs__(x):   return _casadi.fabs(x)
      def __iadd__(x, y):      return _casadi.plus(x, y)
      def __isub__(x, y):      return _casadi.minus(x, y)
      def __imul__(x, y):      return _casadi.times(x, y)
      def __itruediv__(x, y):  return _casadi.rdivide(x, y)
      def __ifloordiv__(x, y): return _casadi.floor(_casadi.rdivide(x, y))
      def __imod__(x, y):      return x - y * _casadi.floor(_casadi.rdivide(x, y))
      def __ipow__(x, n):      return _casadi.power(x, n)
      def __arctan2__(x, y): return _casadi.atan2(x, y)
      def __rarctan2__(y, x): return _casadi.atan2(x, y)
      def fmin(self, other): return _casadi.fmin(self, other)
      def fmax(self, other): return _casadi.fmax(self, other)
      def __fmin__(x, y): return _casadi.fmin(x, y)
      def __rfmin__(y, x): return _casadi.fmin(x, y)
      def __fmax__(x, y): return _casadi.fmax(x, y)
      def __rfmax__(y, x): return _casadi.fmax(x, y)
      def logic_and(self, other): return _casadi.logic_and(self, other)
      def logic_or(self, other): return _casadi.logic_or(self, other)
      def fabs(self): return _casadi.fabs(self)
      def sqrt(self): return _casadi.sqrt(self)
      def sin(self): return _casadi.sin(self)
      def cos(self): return _casadi.cos(self)
      def tan(self): return _casadi.tan(self)
      def arcsin(self): return _casadi.asin(self)
      def arccos(self): return _casadi.acos(self)
      def arctan(self): return _casadi.atan(self)
      def sinh(self): return _casadi.sinh(self)
      def cosh(self): return _casadi.cosh(self)
      def tanh(self): return _casadi.tanh(self)
      def arcsinh(self): return _casadi.asinh(self)
      def arccosh(self): return _casadi.acosh(self)
      def arctanh(self): return _casadi.atanh(self)
      def exp(self): return _casadi.exp(self)
      def log(self): return _casadi.log(self)
      def log10(self): return _casadi.log10(self)
      def log1p(self): return _casadi.log1p(self)
      def expm1(self): return _casadi.expm1(self)
      def floor(self): return _casadi.floor(self)
      def ceil(self): return _casadi.ceil(self)
      def erf(self): return _casadi.erf(self)
      def sign(self): return _casadi.sign(self)
      def fmod(self, other): return _casadi.fmod(self, other)
      def hypot(self, other): return _casadi.hypot(self, other)
      def remainder(self, other): return _casadi.remainder(self, other)
      def __copysign__(x, y): return _casadi.copysign(x, y)
      def __rcopysign__(y, x): return _casadi.copysign(x, y)
      def copysign(self, other): return _casadi.copysign(self, other)
      def rcopysign(self, other): return _casadi.copysign(other, self)
      def __constpow__(x, y): return _casadi.constpow(x, y)
      def __rconstpow__(y, x): return _casadi.constpow(x, y)
      def constpow(self, other): return _casadi.constpow(self, other)
      def rconstpow(self, other): return _casadi.constpow(other, self)
    %}
    /* PEP-484 stubs for the %pythoncode operators above. */
    %stub_CasadiMatrix_binop(__add__)
    %stub_CasadiMatrix_binop(__radd__)
    %stub_CasadiMatrix_binop(__sub__)
    %stub_CasadiMatrix_binop(__rsub__)
    %stub_CasadiMatrix_binop(__mul__)
    %stub_CasadiMatrix_binop(__rmul__)
    %stub_CasadiMatrix_binop(__truediv__)
    %stub_CasadiMatrix_binop(__rtruediv__)
    %stub_CasadiMatrix_binop(__floordiv__)
    %stub_CasadiMatrix_binop(__rfloordiv__)
    %stub_CasadiMatrix_binop(__mod__)
    %stub_CasadiMatrix_binop(__rmod__)
    %stub_overload_method_selftyped(__pow__, DM, DM, other: _DM, modulo: _DM | None = None)
    %stub_overload_method_selftyped(__pow__, SX, SX, other: _SX, modulo: _SX | None = None)
    %stub_overload_method_selftyped(__pow__, MX, MX, other: _MX, modulo: _MX | None = None)
    %stub_CasadiMatrix_binop(__rpow__)
    %stub_CasadiMatrix_binop(__iadd__)
    %stub_CasadiMatrix_binop(__isub__)
    %stub_CasadiMatrix_binop(__imul__)
    %stub_CasadiMatrix_binop(__itruediv__)
    %stub_CasadiMatrix_binop(__ifloordiv__)
    %stub_CasadiMatrix_binop(__imod__)
    %stub_CasadiMatrix_binop(__ipow__)
    %stub_CasadiMatrix_unop(__abs__)
    %stub_CasadiMatrix_unop(__trunc__)
    %stub_CasadiMatrix_unop(__floor__)
    %stub_CasadiMatrix_unop(__ceil__)
    %stub_CasadiMatrix_binop(__lt__)
    %stub_CasadiMatrix_binop(__le__)
    %stub_CasadiMatrix_binop(__gt__)
    %stub_CasadiMatrix_binop(__ge__)
    %stub_CasadiMatrix_binop(__eq__)
    %stub_CasadiMatrix_binop(__ne__)
    %stub_CasadiMatrix_binop(__req__)
    %stub_CasadiMatrix_binop(__rne__)
    %stub_CasadiMatrix_binop(__rlt__)
    %stub_CasadiMatrix_binop(__rle__)
    %stub_CasadiMatrix_binop(__rgt__)
    %stub_CasadiMatrix_binop(__rge__)
    %stub_method(__round__, Self, ndigits: int | None = None)
    %stub_overload_method_selftyped(__divmod__,  %arg(tuple[DM, DM]), DM, other: _DM)
    %stub_overload_method_selftyped(__divmod__,  %arg(tuple[SX, SX]), SX, other: _SX)
    %stub_overload_method_selftyped(__divmod__,  %arg(tuple[MX, MX]), MX, other: _MX)
    %stub_overload_method_selftyped(__rdivmod__, %arg(tuple[DM, DM]), DM, other: _DM)
    %stub_overload_method_selftyped(__rdivmod__, %arg(tuple[SX, SX]), SX, other: _SX)
    %stub_overload_method_selftyped(__rdivmod__, %arg(tuple[MX, MX]), MX, other: _MX)
    /* numpy-style method spellings of the elementwise functions. */
    %stub_CasadiMatrix_unop(fabs)
    %stub_CasadiMatrix_unop(sqrt)
    %stub_CasadiMatrix_unop(sin)
    %stub_CasadiMatrix_unop(cos)
    %stub_CasadiMatrix_unop(tan)
    %stub_CasadiMatrix_unop(arcsin)
    %stub_CasadiMatrix_unop(arccos)
    %stub_CasadiMatrix_unop(arctan)
    %stub_CasadiMatrix_unop(sinh)
    %stub_CasadiMatrix_unop(cosh)
    %stub_CasadiMatrix_unop(tanh)
    %stub_CasadiMatrix_unop(arcsinh)
    %stub_CasadiMatrix_unop(arccosh)
    %stub_CasadiMatrix_unop(arctanh)
    %stub_CasadiMatrix_unop(exp)
    %stub_CasadiMatrix_unop(log)
    %stub_CasadiMatrix_unop(log10)
    %stub_CasadiMatrix_unop(log1p)
    %stub_CasadiMatrix_unop(expm1)
    %stub_CasadiMatrix_unop(floor)
    %stub_CasadiMatrix_unop(ceil)
    %stub_CasadiMatrix_unop(erf)
    %stub_CasadiMatrix_unop(sign)
    %stub_CasadiMatrix_binop(fmin)
    %stub_CasadiMatrix_binop(fmax)
    %stub_CasadiMatrix_binop(logic_and)
    %stub_CasadiMatrix_binop(logic_or)
    %stub_CasadiMatrix_binop(fmod)
    %stub_CasadiMatrix_binop(hypot)
    %stub_CasadiMatrix_binop(remainder)
    %stub_CasadiMatrix_binop(copysign)
    %stub_CasadiMatrix_binop(rcopysign)
    %stub_CasadiMatrix_binop(constpow)
    %stub_CasadiMatrix_binop(rconstpow)
    %stub_CasadiMatrix_binop(__arctan2__)
    %stub_CasadiMatrix_binop(__rarctan2__)
    %stub_CasadiMatrix_binop(__fmin__)
    %stub_CasadiMatrix_binop(__rfmin__)
    %stub_CasadiMatrix_binop(__fmax__)
    %stub_CasadiMatrix_binop(__rfmax__)
    %stub_CasadiMatrix_binop(__copysign__)
    %stub_CasadiMatrix_binop(__rcopysign__)
    %stub_CasadiMatrix_binop(__constpow__)
    %stub_CasadiMatrix_binop(__rconstpow__)
  }

  %extend GenericMatrixCommon {
    %pythoncode %{
      def __mldivide__(x, y): return _casadi.mldivide(x, y)
      def __rmldivide__(y, x): return _casadi.mldivide(x, y)
      def __mrdivide__(x, y): return _casadi.mrdivide(x, y)
      def __rmrdivide__(y, x): return _casadi.mrdivide(x, y)
      def __mpower__(x, y): return _casadi.mpower(x, y)
      def __rmpower__(y, x): return _casadi.mpower(x, y)
    %}
    %stub_CasadiMatrix_binop(__mldivide__)
    %stub_CasadiMatrix_binop(__rmldivide__)
    %stub_CasadiMatrix_binop(__mrdivide__)
    %stub_CasadiMatrix_binop(__rmrdivide__)
    %stub_CasadiMatrix_binop(__mpower__)
    %stub_CasadiMatrix_binop(__rmpower__)
  }

} // namespace casadi

%extend casadi::Callback {
  // eval(arg, /): arg is a sequence of DM
  %stub_method(eval, %arg(Sequence[_DM]), %arg(__arg: Sequence[DM]))
  // eval_buffer(arg, res): tuples of memoryviews
  %stub_method(eval_buffer, int, %arg(__arg: tuple[memoryview, ...], res: tuple[memoryview, ...]))
  %stub_method(get_sparsity_in,  Sparsity, i: int)
  %stub_method(get_sparsity_out, Sparsity, i: int)
}

%extend casadi::CasadiMeta {
  static const char* swig_flags() { return CASADI_SWIG_FLAGS; }
};

%extend casadi::DeserializerBase {
  %stub_method0(unpack, Any)
  %pythoncode %{
    def unpack(self):
      type = SerializerBase.type_to_string(self._pop_type())
      f = getattr(self, "blind_unpack_"+type)
      return f()
  %}
}
%pythoncode %{
veccat.__doc__, vertcat.__doc__, horzcat.__doc__, diagcat.__doc__ = \
    _veccat.__doc__, _vertcat.__doc__, _horzcat.__doc__, _diagcat.__doc__

try:
  import threading
  _thread_local = threading.local()
except:
  threading_available = True
  _thread_local = globals()

def _current_pickle_context():
  return getattr(_thread_local, "casadi_pickle_ctx", None)

def _current_unpickle_context():
  return getattr(_thread_local, "casadi_unpickle_ctx", None)

class global_pickle_context:
    def __enter__(self):
        self.ctx = StringSerializer()
        _thread_local.casadi_pickle_ctx = self.ctx
        return self.ctx

    def __exit__(self, *args):
        _thread_local.casadi_pickle_ctx = None
      
class global_unpickle_context:
    def __enter__(self):
        self.ctx = StringDeserializer("")
        _thread_local.casadi_unpickle_ctx = self.ctx
        return self.ctx

    def __exit__(self, *args):
        _thread_local.casadi_unpickle_ctx = None
%}

// Stubs for pure-pythoncode classes
%stub_class_begin(global_pickle_context)
%stub_method0(__enter__, StringSerializer)
%stub_method(__exit__, None, *args: object)

%stub_class_begin(global_unpickle_context)
%stub_method0(__enter__, StringDeserializer)
%stub_method(__exit__, None, *args: object)


%extend casadi::Opti {
  %pythoncode %{

    @staticmethod
    def _callback(self,fh=None):
      if fh is None:
        self.callback_class();
        return
      class OptiCallbackHelper(OptiCallback):
          def __init__(self, callback):
            OptiCallback.__init__(self)
            self.callback = callback

          def call(self, i):
            self.callback(i)

      self._fh = fh
      self._cb = OptiCallbackHelper(fh);
      self.callback_class(self._cb);


    def callback(self,fh=None):
      self._callback(self,fh)


  %}
  // Opti.callback (%pythoncode); fh gets the iteration index
  %stub_method(callback, None, %arg(fh: Callable[[int], Any] | None = ...))

}

namespace casadi{
%extend PrintableCommon {
  %pythoncode %{
    def __str__(self): return self.str()
    def repr(self): return self.type_name() + '(' + self.str() + ')'
    # Value-copy semantics for all casadi types (DM/SX/MX/Sparsity/Function/
    # ...): copy via the copy constructor.  This base injects __copy__ /
    # __deepcopy__ into every subclass, replacing the old global-`object`
    # shadow hack.
    def __copy__(self): return self.__class__(self)
    def __deepcopy__(self, memo=None): return self.__class__(self)
  %}
  %stub_method0(repr, builtins.str)
  %stub_method0(__copy__, Self)
  %stub_method(__deepcopy__, Self, memo: dict[int, Any] | None = None)
}
} // namespace casadi

%extend casadi::MX{
  %python_array_wrappers(1002.0, Any)
}

namespace casadi {
%extend Matrix<SXElem>{
  %python_array_wrappers(1001.0, Any)
}
} // namespace casadi
