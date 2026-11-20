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

// Host conversions: is_null, full/sparse, casadi_host_conv overloads
%fragment("casadi_host", "header", fragment="casadi_decl", fragment="JuliaArrayTraits") {
  namespace casadi {
    bool is_null(GUESTOBJECT *p) {
      return !p || jl_is_nothing(p);
    }

    GUESTOBJECT* full(const DM& m, bool simplify) {
      return 0;
    }

    // Convert to a sparse matrix
    GUESTOBJECT* sparse(const DM& m) {
      return 0;
    }
  } // namespace casadi

  namespace swig {
    // asval through in_place: container elements and map values convert into their slot
    template<class M> struct casadi_asval {
      static int asval(GUESTOBJECT* p, M* v) {
        if (!v) return traits_asptr<M>::asptr(p, 0);
        M* vp = v;
        if (!traits_asptr<M>::in_place(p, &vp)) return SWIG_ERROR;
        if (vp != v) *v = *vp;
        return SWIG_OK;
      }
    };

      // Real Julia array of at most two dimensions to dense DM (vector -> column)
    inline bool from_julia_array(jl_value_t* p, casadi::DM* m) {
        using namespace casadi;
        jl_value_t* aa = jl_get_global(jl_base_module, jl_symbol("AbstractArray"));
        if (!aa || !jl_isa(p, aa)) return false;
        jl_value_t *et = 0, *sz = 0, *arr = 0, *atype = 0;
        JL_GC_PUSH4(&et, &sz, &arr, &atype);
        et = jl_call1(jl_get_function(jl_base_module, "eltype"), p);
        sz = jl_call1(jl_get_function(jl_base_module, "size"), p);
        size_t nd = sz && jl_is_tuple(sz) ? jl_nfields(sz) : 3;
        bool ok = et && nd <= 2;
        casadi_int dims[2] = {1, 1};
        for (size_t i = 0; ok && i < nd; ++i) dims[i] = jl_unbox_int64(jl_fieldref(sz, i));
        bool real = ok && jl_is_type(et)
          && jl_subtype(et, jl_get_global(jl_base_module, jl_symbol("Real")));
        if (ok && !real) {
          // `[]` has eltype Any: accept it as an empty DM, reject anything else
          ok = dims[0] * dims[1] == 0;
          if (ok && m) *m = DM(0, 0);
        } else if (ok && m) {
          atype = jl_apply_array_type((jl_value_t*)jl_float64_type, nd);
          arr = jl_call2(jl_get_function(jl_base_module, "convert"), atype, p);
          jl_value_t* ptr = arr ? jl_call1(jl_get_function(jl_base_module, "pointer"), arr) : 0;
          ok = ptr != 0;
          if (ok) {
            const double* d = static_cast<const double*>(jl_unbox_voidpointer(ptr));
            *m = DM::zeros(dims[0], dims[1]);
            std::copy(d, d + dims[0] * dims[1], m->ptr());
          }
        }
        JL_GC_POP();
        return ok;
      }

    inline int casadi_host_conv(GUESTOBJECT* p, casadi::DM* m) {
      if (from_julia_array(p, m)) return true;

      // No match
      return false;
    }

    // Returns as the Julia value, not the proxy
    template<> struct julia_eltype<casadi::GenericType> {
      static jl_value_t* get() { return (jl_value_t*)jl_any_type; }
    };

    inline GUESTOBJECT* casadi_host_from(const casadi::GenericType& a) {
      using namespace casadi;
      switch (a.getType()) {
      case OT_NULL:
      case OT_VOIDPTR:
        return jl_nothing;
      default: return 0;
      }
    }

    // asval through in_place for every wrapped type
    template<> struct traits_asval<casadi::Sparsity> : casadi_asval<casadi::Sparsity> {};
    template<> struct traits_asval<casadi::Function> : casadi_asval<casadi::Function> {};
    template<> struct traits_asval<casadi::Slice> : casadi_asval<casadi::Slice> {};
    template<> struct traits_asval<casadi::Matrix<double> > : casadi_asval<casadi::Matrix<double> > {};
    template<> struct traits_asval<casadi::Matrix<casadi_int> > : casadi_asval<casadi::Matrix<casadi_int> > {};
    template<> struct traits_asval<casadi::Matrix<casadi::SXElem> > : casadi_asval<casadi::Matrix<casadi::SXElem> > {};
    template<> struct traits_asval<casadi::SXElem> : casadi_asval<casadi::SXElem> {};
    template<> struct traits_asval<casadi::MX> : casadi_asval<casadi::MX> {};
    template<> struct traits_asval<casadi::GenericType> : casadi_asval<casadi::GenericType> {};
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
