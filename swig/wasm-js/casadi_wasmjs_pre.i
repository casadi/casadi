// Docstrings keep the \lang groups of the target language
%feature("docstring:lang", "javascript");

// Director typemaps handing wrapped proxies to JS

%define %casadi_director_class(xType...)
// SWIG_NewPointerObj wraps via clientdata
%typemap(directorin, noblock=1) xType, const xType, const xType& {
  {
    xType *__p = new xType($1);
    $input = SWIG_NewPointerObj(__p, $descriptor(xType *), 0);
  }
}
%typemap(directorout, noblock=1) xType, const xType {
  {
    emscripten::val __r = emscripten::val::take_ownership($input);
    EM_VAL __h = __r.release_ownership();
    void *__p = 0;
    (void)SWIG_WASMJS_ConvertPtr(__h, &__p, 0, 0);
    $result = __p ? *static_cast< xType *>(__p) : xType();
  }
}
%enddef

%define %casadi_director_vec(xElemType...)
%typemap(directorin, noblock=1) std::vector< xElemType >,
                                 const std::vector< xElemType >&,
                                       std::vector< xElemType >& {
  // SWIG_NewPointerObj wraps via clientdata
  {
    emscripten::val __a = emscripten::val::array();
    for (size_t __i = 0; __i < $1.size(); ++__i) {
      xElemType *__ep = new xElemType($1[__i]);
      EM_VAL __eh = SWIG_NewPointerObj(__ep, $descriptor(xElemType *), 0);
      __a.call<void>("push", emscripten::val::take_ownership(__eh));
    }
    $input = __a.release_ownership();
  }
}
%typemap(directorout, noblock=1) std::vector< xElemType > {
  {
    emscripten::val __r = emscripten::val::take_ownership($input);
    std::vector< xElemType > __v;
    size_t __n = __r["length"].as<size_t>();
    for (size_t __i = 0; __i < __n; ++__i) {
      emscripten::val __it = __r[__i];
      EM_VAL __h = __it.release_ownership();
      void *__p = 0;
      (void)SWIG_WASMJS_ConvertPtr(__h, &__p, 0, 0);
      if (__p) __v.push_back(*static_cast< xElemType *>(__p));
    }
    $result = __v;
  }
}
%enddef

// std::vector<std::string>: elements are plain JS strings
%define %casadi_director_vec_str
%typemap(directorin, noblock=1) std::vector< std::string >,
                                 const std::vector< std::string >&,
                                       std::vector< std::string >& {
  {
    emscripten::val __a = emscripten::val::array();
    for (size_t __i = 0; __i < $1.size(); ++__i) {
      __a.call<void>("push", emscripten::val(std::string($1[__i])));
    }
    $input = __a.release_ownership();
  }
}
%typemap(directorout, noblock=1) std::vector< std::string > {
  {
    emscripten::val __r = emscripten::val::take_ownership($input);
    std::vector< std::string > __v;
    size_t __n = __r["length"].as<size_t>();
    for (size_t __i = 0; __i < __n; ++__i) {
      __v.push_back(__r[__i].as<std::string>());
    }
    $result = __v;
  }
}
%enddef

// jsout for std::map<std::string, X>: unbox into a native JS object
%define %casadi_dict_jsout(xValueType...)
%typemap(jsout) std::map< std::string, xValueType >,
                 const std::map< std::string, xValueType >&,
                       std::map< std::string, xValueType >&
"M.__swig_release_handle($call)"
%enddef

%fragment("casadi_containers", "header", fragment="casadi_aux", fragment="StdVectorTraits",
          fragment="casadi_map_traits", fragment="StdPairTraits", fragment=SWIG_Traits_frag(bool),
          fragment=SWIG_Traits_frag(long long), fragment=SWIG_Traits_frag(double),
          fragment=SWIG_Traits_frag(std::string)) {
  namespace swig {
    // Numbers also accept a JS boolean
    template<class T> struct casadi_bool_number {
      static bool asref(EM_VAL p, T** m) {
        if (SWIG_IsOK(asval(p, *m))) return true;
        bool b;
        if (!SWIG_IsOK(asval(p, &b))) return false;
        **m = b;
        return true;
      }
      static bool check(EM_VAL p) {
        return SWIG_IsOK(asval(p, (T*)0)) || SWIG_IsOK(asval(p, (bool*)0));
      }
    };
    template<> struct traits_asref<double, value_category, void> : casadi_bool_number<double> {};
    template<> struct traits_check<double, value_category> : casadi_bool_number<double> {};
    template<> struct traits_asref<long long, value_category, void> : casadi_bool_number<long long> {};
    template<> struct traits_check<long long, value_category> : casadi_bool_number<long long> {};

    // A null or undefined dict is an empty dict
    template<class T, class C, class A> struct traits_asref<std::map<std::string, T, C, A>, pointer_category, void> {
      typedef std::map<std::string, T, C, A> map_type;
      static bool asref(EM_VAL p, map_type** m) {
        if (casadi::is_null(p)) return true;
        map_type* r = 0;
        int res = asptr(p, &r);
        if (!SWIG_IsOK(res) || !r) return false;
        if (SWIG_IsNewObj(res)) {
          **m = std::move(*r);
          delete r;
        } else {
          *m = r;
        }
        return true;
      }
    };
    template<class T, class C, class A> struct traits_check<std::map<std::string, T, C, A>, pointer_category> {
      static bool check(EM_VAL p) {
        return casadi::is_null(p) || SWIG_IsOK(asptr(p, (std::map<std::string, T, C, A>**)0));
      }
    };

    // Array elements convert with the casadi rules; in a check, a wrapped element must be of the element type itself
    template<class T, class A> struct JSArrayProtocol<std::vector<T, A>, T> {
      static bool assign(const emscripten::val& arr, std::vector<T, A>* seq) {
        unsigned n = arr["length"].as<unsigned>();
        seq->reserve(n);
        T tmp;
        for (unsigned i = 0; i < n; ++i) {
          emscripten::val el = arr[i];
          T* e = &tmp;
          if (!traits_asref<T>::asref(el.as_handle(), &e)) return false;
          seq->push_back(*e);
        }
        return true;
      }
      static bool exact(EM_VAL p, pointer_category) {
        T* r = 0;
        int res = asptr(p, &r);
        if (SWIG_IsOK(res) && SWIG_IsNewObj(res)) delete r;
        return SWIG_IsOK(res) && !SWIG_IsNewObj(res);
      }
      static bool exact(EM_VAL, value_category) { return false; }
      static bool check(const emscripten::val& arr) {
        unsigned n = arr["length"].as<unsigned>();
        for (unsigned i = 0; i < n; ++i) {
          emscripten::val el = arr[i];
          std::string ty = el.typeOf().as<std::string>();
          bool tagged = (ty == "object" || ty == "function") && !el.isNull()
            && el["_swig_type"].typeOf().as<std::string>() == "string";
          if (tagged) {
            if (!exact(el.as_handle(), typename traits<T>::category())) return false;
          } else {
            T tmp, *e = &tmp;
            if (!swig::asref(el.as_handle(), &e)) return false;
          }
        }
        return true;
      }
    };
  } // namespace swig
}

// Host conversions: is_null, full/sparse, casadi_host_conv overloads
%fragment("casadi_host", "header", fragment="casadi_decl") {
  namespace casadi {
    bool is_null(GUESTOBJECT *p) {
      if (!p) return true;
      emscripten::val v = SWIG_WASMJS_Borrow(p);
      return v.isNull() || v.isUndefined();
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
      // A JS number or BigInt as a double
    inline bool js_number(const emscripten::val& el, double* d) {
        std::string ty = el.typeOf().as<std::string>();
        if (ty == "number") *d = el.as<double>();
        else if (ty == "bigint") *d = emscripten::val::global("Number")(el).as<double>();
        else return false;
        return true;
      }

    inline int casadi_host_conv(GUESTOBJECT* p, casadi::DM* m) {
      // JS Array of numbers -> column vector; nested [[a,b],[c,d]] -> dense matrix, row-major in JS
      emscripten::val v = SWIG_WASMJS_Borrow(p);
      if (!v.isArray()) return false;
      unsigned n = v["length"].as<unsigned>();
      bool nested = false;
      if (n > 0) {
        emscripten::val first = v[0];
        nested = first.isArray();
      }
      if (nested) {
        unsigned ncols = 0;
        for (unsigned i = 0; i < n; ++i) {
          emscripten::val row = v[i];
          if (!row.isArray()) return false;
          unsigned w = row["length"].as<unsigned>();
          if (i == 0) ncols = w;
          else if (w != ncols) return false;
        }
        std::vector<double> data((size_t)n * (size_t)ncols);
        for (unsigned i = 0; i < n; ++i) {
          emscripten::val row = v[i];
          for (unsigned j = 0; j < ncols; ++j) {
            emscripten::val el = row[j];
            // Column-major: data[j*n + i]
            if (!js_number(el, &data[(size_t)j * (size_t)n + (size_t)i])) return false;
          }
        }
        if (m) *m = casadi::Matrix<double>::reshape(casadi::Matrix<double>(data), n, ncols);
        return true;
      }
      std::vector<double> data(n);
      for (unsigned i = 0; i < n; ++i) {
        emscripten::val el = v[i];
        if (!js_number(el, &data[i])) return false;
      }
      if (m) *m = casadi::Matrix<double>(data);
      return true;
    }
  } // namespace swig
}

// Target typemaps of each casadi-typemapped type
%define CASADI_TARGET_INPUT_TYPEMAPS(xName, xStubIn, xTsStub, xType...)
 // Wrapped types cross the wasm boundary as EM_VAL
%typemap(ctype) xType, const xType&, xType& "EM_VAL"
 // Non-const reference: standard SWIG convention
%typemap(in, doc=xName, pystub_in=xStubIn, tsstub_in=xTsStub, noblock=1) xType & {
  if (!SWIG_IsOK(SWIG_ConvertPtr($input, (void**)&$1, $1_descriptor, 0))) {
    SWIG_exception_fail(SWIG_TypeError, "Failed to convert input $argnum to type '" xName "'.");
  }
 }
%typemap(freearg, noblock=1) xType & {}
%enddef

// Named std::vector class
%define CASADI_TARGET_TEMPLATE_VEC(xJSName, xElemType...)
%template(xJSName) std::vector< xElemType >;
%wasm_vec(xElemType, xJSName)
%enddef
