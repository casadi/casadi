// Docstrings keep the \lang groups of the target language
%feature("docstring:lang", "julia");

%{
  namespace casadi {
    // Printout goes to Julia's stdout/stderr
    static jl_value_t* julia_write = 0;
    static jl_value_t* julia_flush = 0;
    static void julialogger(const char* s, std::streamsize num, bool error) {
      // Threads Julia does not know print directly
      if (!jl_get_pgcstack()) return casadi::Logger::writeDefault(s, num, error);
      SwigJlGcUnsafe enter;
      jl_value_t* str = jl_pchar_to_string(s, num);
      JL_GC_PUSH1(&str);
      if (!jl_call2(julia_write, str, error ? jl_true : jl_false)) jl_exception_clear();
      JL_GC_POP();
    }
    static void juliaflush(bool error) {
      if (!jl_get_pgcstack()) return casadi::Logger::flushDefault(error);
      SwigJlGcUnsafe enter;
      if (!jl_call1(julia_flush, error ? jl_true : jl_false)) jl_exception_clear();
    }
  } // namespace casadi
%}
%wrapper %{
// Set from the Julia module's __init__; the module keeps write and flush alive
extern "C" SWIGEXPORT void casadi_julia_logger(jl_value_t* write, jl_value_t* flush) {
  casadi::julia_write = write;
  casadi::julia_flush = flush;
  casadi::Logger::writeFun = casadi::julialogger;
  casadi::Logger::flush = casadi::juliaflush;
}
%}

%fragment("casadi_containers", "header", fragment="casadi_aux", fragment="StdVectorTraits",
          fragment="casadi_map_traits", fragment="StdPairTraits", fragment=SWIG_Traits_frag(bool),
          fragment=SWIG_Traits_frag(long long), fragment=SWIG_Traits_frag(double),
          fragment=SWIG_Traits_frag(std::string)) {
  namespace swig {
    // Other AbstractVectors (ranges, views, ...): collect, then convert the Array
    template<class T> struct traits_asptr<std::vector<T, std::allocator<T> > > {
      static int asptr(jl_value_t* p, std::vector<T>** m) {
        jl_value_t* av = p && !jl_is_array(p) ? jl_get_global(jl_base_module, jl_symbol("AbstractVector")) : 0;
        if (!av || !jl_isa(p, av)) return traits_asptr_stdseq<std::vector<T> >::asptr(p, m);
        jl_value_t* c = 0;
        JL_GC_PUSH1(&c);
        c = jl_call1(jl_get_function(jl_base_module, "collect"), p);
        int res = c && jl_is_array(c) ? traits_asptr_stdseq<std::vector<T> >::asptr(c, m) : SWIG_ERROR;
        JL_GC_POP();
        return res;
      }
    };
  } // namespace swig
}

// Target typemaps of each casadi-typemapped type
%define CASADI_TARGET_INPUT_TYPEMAPS(xName, xStubIn, xTsStub, xType...)
%julia_value_type(xType)
%enddef

// Returns follow the Lib convention
%define CASADI_TARGET_OUTPUT_TYPEMAPS(xName, xStubIn, xStubOut, xTsStub, xType...)
%enddef

// Type names in error messages and docs
  #define L_INT "Int"
  #define L_BOOL "Bool"
  #define LPAIR(A,B) "Tuple{" A "," B "}"
  #define L_DOUBLE "Float64"
  #define L_DICT "Dict"
  #define LDICT(ARG) "Dict{String," ARG "}"
  #define LL "Vector{"
  #define LR "}"
  #define L_STR "String"
  #define LABEL LPAIR(L_INT, L_INT)

%define make_property(class, name, type)
// Julia: obj.name is the value itself (opti.f, opti.x, sol.debug), not a bound method
  %feature("julia:property") class ## :: ## name;
%enddef
