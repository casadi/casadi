// Overrides of stock fragments; fragments are first-registration-wins, so this precedes stl.i

// Descriptor of p if p is a DM, SX or MX but not a subclass, else 0; these define no conversion hooks
%fragment("casadi_exact_type", "header") {
namespace casadi {
  SWIGINTERN swig_type_info* exact_type(PyObject* p) {
    static swig_type_info* const types[] = {$descriptor(casadi::Matrix<double>*),
      $descriptor(casadi::Matrix<casadi::SXElem>*), $descriptor(casadi::MX*)};
    for (swig_type_info* t : types) {
      SwigPyClientData* cd = static_cast<SwigPyClientData*>(t->clientdata);
      if (cd && reinterpret_cast<PyObject*>(Py_TYPE(p)) == cd->klass) return t;
    }
    return 0;
  }
} // namespace casadi
}

// New reference to p.item() for a size-1 numpy-like object, else 0 (no error left set)
%fragment("casadi_np_item", "header", fragment="casadi_exact_type") {
namespace casadi {
  SWIGINTERN PyObject* np_item(PyObject* p) {
    if (PyErr_Occurred()) PyErr_Clear();
    if (exact_type(p)) return 0;
    if (!PyObject_HasAttrString(p, "__array__")) return 0;
    PyObject* sz = PyObject_GetAttrString(p, "size");
    if (!sz) { PyErr_Clear(); return 0; }
    long long n = -1;
    int res = SWIG_AsVal(long long)(sz, &n);
    Py_DECREF(sz);
    if (!SWIG_IsOK(res) || n != 1) return 0;
    PyObject* it = PyObject_CallMethod(p, (char*) "item", 0);
    if (!it) PyErr_Clear();
    return it;
  }
} // namespace casadi
}

%fragment(SWIG_Traits_frag(double), "header", fragment=SWIG_AsVal_frag(double),
          fragment=SWIG_From_frag(double), fragment=SWIG_AsVal_frag(long long),
          fragment="StdTraits", fragment="casadi_np_item") {
namespace swig {
  template <> struct traits< double > {
    typedef value_category category;
    static const char* type_name() { return "double"; }
  };
  template <> struct traits_asval< double > {
    static int asval(PyObject* obj, double* val) {
      int res = SWIG_AsVal(double)(obj, val);
      if (SWIG_IsOK(res)) return res;
      PyObject* it = casadi::np_item(obj);
      if (!it) return SWIG_ERROR;
      res = asval(it, val);
      Py_DECREF(it);
      return res;
    }
  };
  template <> struct traits_from< double > {
    static PyObject* from(const double& val) { return SWIG_From(double)(val); }
  };
}
}

%fragment(SWIG_Traits_frag(long long), "header", fragment=SWIG_AsVal_frag(long long),
          fragment=SWIG_From_frag(long long), fragment="StdTraits", fragment="casadi_np_item") {
namespace swig {
  template <> struct traits< long long > {
    typedef value_category category;
    static const char* type_name() { return "long long"; }
  };
  template <> struct traits_asval< long long > {
    static int asval(PyObject* obj, long long* val) {
      int res = SWIG_AsVal(long long)(obj, val);
      if (SWIG_IsOK(res)) return res;
      PyObject* it = casadi::np_item(obj);
      if (!it) return SWIG_ERROR;
      res = asval(it, val);
      Py_DECREF(it);
      return res;
    }
  };
  template <> struct traits_from< long long > {
    static PyObject* from(const long long& val) { return SWIG_From(long long)(val); }
  };
}
}

// std::pair: a 2-tuple only (stock also takes a 2-element sequence)
%fragment("StdPairTraits", "header", fragment="StdTraits") {
namespace swig {
  template <class T, class U> struct traits_asptr<std::pair<T, U> > {
    static int asptr(PyObject* obj, std::pair<T, U>** val) {
      if (!PyTuple_Check(obj) || PyTuple_Size(obj) != 2) return SWIG_ERROR;
      PyObject* first = PyTuple_GetItem(obj, 0);
      PyObject* second = PyTuple_GetItem(obj, 1);
      if (!val) {
        if (SWIG_IsOK(swig::asval(first, (T*)0)) && SWIG_IsOK(swig::asval(second, (U*)0))) return SWIG_OK;
        if (PyErr_Occurred()) PyErr_Clear();
        return SWIG_ERROR;
      }
      std::pair<T, U>* p = new std::pair<T, U>();
      if (!SWIG_IsOK(swig::asval(first, &p->first)) || !SWIG_IsOK(swig::asval(second, &p->second))) {
        delete p;
        if (PyErr_Occurred()) PyErr_Clear();
        return SWIG_ERROR;
      }
      *val = p;
      return SWIG_NEWOBJ;
    }
  };
  template <class T, class U> struct traits_from<std::pair<T, U> > {
    static PyObject* from(const std::pair<T, U>& val) {
      PyObject* first = swig::from(val.first);
      if (!first) return 0;
      PyObject* second = swig::from(val.second);
      if (!second) { Py_DECREF(first); return 0; }
      PyObject* obj = PyTuple_New(2);
      PyTuple_SetItem(obj, 0, first);
      PyTuple_SetItem(obj, 1, second);
      return obj;
    }
  };
}
}


#ifdef CASADI_WITH_PYTHON_GIL_RELEASE
%{
  // This .cxx was swig-compiled with WITH_PYTHON_GIL_RELEASE option
  #define CASADI_WITH_PYTHON_GIL_RELEASE
%}
#else //CASADI_WITH_PYTHON_GIL_RELEASE
%{
  // This .cxx was swig-compiled without WITH_PYTHON_GIL_RELEASE option
  #undef CASADI_WITH_PYTHON_GIL_RELEASE
%}
#endif //CASADI_WITH_PYTHON_GIL_RELEASE

%ignore CASADI_SWIG_FLAGS;
%include "swig_config.h"

%{
  namespace casadi {

    // Redirect printout
    static void pythonlogger(const char* s, std::streamsize num, bool error) {
#ifndef CASADI_WITH_PYTHON_GIL_RELEASE
      if (!casadi::InterruptHandler::is_main_thread()) {
        casadi::Logger::writeDefault(s, num, error);
        return;
      }
#endif // CASADI_WITH_PYTHON_GIL_RELEASE
      int n = num;
#ifdef CASADI_WITH_PYTHON_GIL_RELEASE
      SWIG_PYTHON_THREAD_BEGIN_BLOCK;
#endif // CASADI_WITH_PYTHON_GIL_RELEASE
      while (n>0) {
        if (error) {
          PySys_WriteStderr("%.*s", std::min(n, 1000), s);
        } else {
          PySys_WriteStdout("%.*s", std::min(n, 1000), s);
        }
        n -= 1000;
        s += 1000;
      }
#ifdef CASADI_WITH_PYTHON_GIL_RELEASE
      SWIG_PYTHON_THREAD_END_BLOCK;
#endif // CASADI_WITH_PYTHON_GIL_RELEASE
    }

    static bool pythoncheckinterrupted() {
      if (!casadi::InterruptHandler::is_main_thread()) return false;
#ifdef CASADI_WITH_PYTHON_GIL_RELEASE
      SWIG_PYTHON_THREAD_BEGIN_BLOCK;
#endif // CASADI_WITH_PYTHON_GIL_RELEASE
      return PyErr_CheckSignals();
      // SWIG_PYTHON_THREAD_END_BLOCK is not needed, destructor will release GIL
    }

    // Decodes to "" rather than dereferencing NULL when str_py is absent or not a str
    std::string python_string_to_std_string(PyObject *str_py) {
      if (!str_py) return "";
#if SWIG_VERSION < 0x040200
      const char *chars = SWIG_Python_str_AsChar(str_py);
      std::string str(chars ? chars : "");
      SWIG_Python_str_DelForPy3(chars);
#else
      PyObject *bytes = NULL;
      const char *chars = SWIG_PyUnicode_AsUTF8AndSize(str_py, NULL, &bytes);
      std::string str(chars ? chars : "");
      Py_XDECREF(bytes);
#endif
      if (!chars) PyErr_Clear();
      return str;
    }

    void handle_director_exception() {
	    std::string msg = "Exception in SWIG director ";
      // Note: CASADI_WITH_PYTHON_GIL_RELEASE case has SWIG_PYTHON_THREAD_BEGIN_BLOCK in the caller
#ifndef CASADI_WITH_PYTHON_GIL_RELEASE
      SWIG_PYTHON_THREAD_BEGIN_BLOCK;
#endif // CASADI_WITH_PYTHON_GIL_RELEASE
      if (PyErr_ExceptionMatches(PyExc_KeyboardInterrupt)) {
        PyErr_Clear();
#ifndef CASADI_WITH_PYTHON_GIL_RELEASE
        SWIG_PYTHON_THREAD_END_BLOCK;
#endif // CASADI_WITH_PYTHON_GIL_RELEASE
        throw casadi::KeyboardInterruptException();
      }
      PyObject *ptype, *pvalue, *ptraceback;
      PyErr_Fetch(&ptype, &pvalue, &ptraceback);
      PyObject* msg_py = PyObject_Str(pvalue);
      msg = python_string_to_std_string(msg_py);
      Py_XDECREF(msg_py);
      PyErr_Restore(ptype, pvalue, ptraceback);
      PyErr_Print();
#ifndef CASADI_WITH_PYTHON_GIL_RELEASE
      SWIG_PYTHON_THREAD_END_BLOCK;
#endif // CASADI_WITH_PYTHON_GIL_RELEASE
      casadi_error(msg.c_str());
	  }
  }

%}
%init %{
  // Set logger functions
  casadi::Logger::writeFun = casadi::pythonlogger;

  // @jgillis: please document
  casadi::InterruptHandler::checkInterrupted = casadi::pythoncheckinterrupted;

  casadi::InterruptHandler::is_main_thread();

%}

// Collect multiple outputs in a tuple (swig/swig#1933)
#undef SWIG_AppendOutput
#if SWIG_VERSION >= 0x040300
#define SWIG_AppendOutput(result, obj)  SWIG_Python_AppendOutput_Tuple(result, obj, $isvoid)
#else
// $isvoid is SWIG >= 4.3; 1 reproduces the older SWIG_Python_AppendOutput
#define SWIG_AppendOutput(result, obj)  SWIG_Python_AppendOutput_Tuple(result, obj, 1)
#endif

%{
SWIGINTERN PyObject*
SWIG_Python_AppendOutput_Tuple(PyObject* result, PyObject* obj, int is_void) {
  PyObject*   o2;
  PyObject*   o3;
  if (!result) {
    result = obj;
  } else if (result == Py_None && is_void) {
    Py_DECREF(result);
    result = obj;
  } else {
    if (!PyTuple_Check(result)) {
      o2 = result;
      result = PyTuple_New(1);
      PyTuple_SET_ITEM(result, 0, o2);
    }
    o3 = PyTuple_New(1);
    PyTuple_SET_ITEM(o3, 0, obj);
    o2 = result;
    result = PySequence_Concat(o2, o3);
    Py_DECREF(o2);
    Py_DECREF(o3);
  }
  return result;
}
%}

%fragment("casadi_containers", "header", fragment="casadi_aux", fragment="StdVectorTraits",
          fragment="casadi_map_traits", fragment="StdPairTraits", fragment=SWIG_Traits_frag(bool),
          fragment=SWIG_Traits_frag(long long), fragment=SWIG_Traits_frag(double),
          fragment=SWIG_Traits_frag(std::string), fragment="casadi_exact_type") {
  namespace swig {
    template<class T, class A> struct traits<std::vector<T, A> > {
      typedef pointer_category category;
      static const char* type_name() {
        static std::string n = "std::vector<" + std::string(swig::type_name<T>()) + " >";
        return n.c_str();
      }
    };
    template<class K, class T, class C, class A> struct traits<std::map<K, T, C, A> > {
      typedef pointer_category category;
      static const char* type_name() {
        static std::string n = "std::map<" + std::string(swig::type_name<K>()) + ","
          + std::string(swig::type_name<T>()) + " >";
        return n.c_str();
      }
    };
    template<class T, class U> struct traits<std::pair<T, U> > {
      typedef pointer_category category;
      static const char* type_name() {
        static std::string n = "std::pair<" + std::string(swig::type_name<T>()) + ","
          + std::string(swig::type_name<U>()) + " >";
        return n.c_str();
      }
    };

    // Not a sequence: str/bytes/dict/set, objects with a conversion hook, non-1-D shape
    inline bool casadi_not_a_sequence(PyObject* p) {
      if (PyList_CheckExact(p) || PyTuple_CheckExact(p)) return false;
      if (casadi::exact_type(p)) return true;
      if (PyDict_Check(p) || PyBytes_Check(p) || PySet_Check(p) || PyUnicode_Check(p)) return true;
      if (PyObject_HasAttrString(p, "__MX__") || PyObject_HasAttrString(p, "__SX__")
          || PyObject_HasAttrString(p, "__DM__")) return true;
      if (PyErr_Occurred()) PyErr_Clear();
      if (PyObject_HasAttrString(p, "shape")) {
        PyObject* shape = PyObject_GetAttrString(p, "shape");
        if (!shape) { PyErr_Clear(); return true; }
        bool is_1d = PyTuple_Check(shape) && PyTuple_Size(shape)==1;
        Py_DECREF(shape);
        if (!is_1d) return true;
      }
      return false;
    }

    template<class T, class A> struct IteratorProtocol<std::vector<T, A>, T> {
      static void assign(PyObject* p, std::vector<T, A>* seq) {
        if (casadi_not_a_sequence(p)) throw std::invalid_argument("not a sequence");
        Py_ssize_t n = PyObject_Length(p);
        if (n > 0) seq->reserve(n);
        if (n < 0) PyErr_Clear();
        SwigVar_PyObject it = PyObject_GetIter(p);
        if (!it) return;
        T tmp;
        SwigVar_PyObject item = PyIter_Next(it);
        while (item) {
          T* e = &tmp;
          if (!traits_asref<T>::asref(item, &e)) throw std::invalid_argument("bad type");
          seq->push_back(*e);
          item = PyIter_Next(it);
        }
      }
      static bool check(PyObject* p) {
        if (casadi_not_a_sequence(p)) return false;
        SwigVar_PyObject it = PyObject_GetIter(p);
        if (!it) return casadi_ok(false);
        SwigVar_PyObject item = PyIter_Next(it);
        while (item) {
          if (!swig::check<T>(item)) return false;
          item = PyIter_Next(it);
        }
        return casadi_ok(!PyErr_Occurred());
      }
    };

    // A sequence converts into the caller's vector in place
    template<class T, class A> struct traits_asref<std::vector<T, A>, pointer_category, void> {
      static bool asref(PyObject* p, std::vector<T, A>** m) {
        if (p == Py_None || SWIG_Python_GetSwigThis(p)) return false;
        try {
          (*m)->clear();
          IteratorProtocol<std::vector<T, A>, T>::assign(p, *m);
        } catch (std::exception&) {
          return casadi_ok(false);
        }
        return casadi_ok(!PyErr_Occurred());
      }
    };

    // std::vector maps to a Python list
    template<class T, class A> struct traits_from_stdseq<std::vector<T, A>, T> {
      static PyObject* from(const std::vector<T, A>& seq) {
        PyObject* ret = PyList_New(seq.size());
        if (!ret) return 0;
        for (size_t k=0; k<seq.size(); ++k) {
          PyObject* el = swig::from<T>(seq[k]);
          if (!el) {
            Py_DECREF(ret);
            return 0;
          }
          PyList_SetItem(ret, k, el);
        }
        return ret;
      }
    };

    // std::map<std::string, T> maps to a dict with str keys
    template<class T, class C, class A> struct traits_asptr<std::map<std::string, T, C, A> > {
      typedef std::map<std::string, T, C, A> map_type;
      // m==0: check only
      static bool fill(PyObject* p, map_type* m) {
        if (!PyDict_Check(p)) return false;
        if (m) m->clear();
        PyObject *key, *value;
        Py_ssize_t pos = 0;
        while (PyDict_Next(p, &pos, &key, &value)) {
          if (!PyUnicode_Check(key)) return false;
          if (m) {
            T& v = (*m)[casadi::python_string_to_std_string(key)];
            T* vp = &v;
            if (!traits_asref<T>::asref(value, &vp)) return false;
            if (vp != &v) v = *vp;
          } else if (!swig::check<T>(value)) {
            return false;
          }
        }
        return true;
      }
      static int asptr(PyObject* p, map_type** m) {
        if (!m) return fill(p, 0) ? SWIG_OK : SWIG_ERROR;
        map_type* r = new map_type();
        if (!fill(p, r)) {
          delete r;
          return SWIG_ERROR;
        }
        *m = r;
        return SWIG_NEWOBJ;
      }
      static bool in_place(PyObject* p, map_type** m) {
        return fill(p, *m);
      }
    };
    template<class T, class C, class A> struct traits_from<std::map<std::string, T, C, A> > {
      static PyObject* from(const std::map<std::string, T, C, A>& a) {
        PyObject *p = PyDict_New();
        if (!p) return 0;
        for (typename std::map<std::string, T, C, A>::const_iterator it=a.begin(); it!=a.end(); ++it) {
          PyObject * e = swig::from(it->second);
          if (!e) {
            Py_DECREF(p);
            return 0;
          }
          PyDict_SetItemString(p, it->first.c_str(), e);
          Py_DECREF(e);
        }
        return p;
      }
    };
  } // namespace swig
}

// Host conversions: is_null, full/sparse, casadi_host_conv overloads
%fragment("casadi_map_traits", "header") {} // Python converts std::map in casadi_containers

%fragment("casadi_host", "header", fragment="casadi_decl", fragment="casadi_exact_type",
          fragment=SWIG_AsVal_frag(double)) {
  namespace casadi {
    bool is_null(GUESTOBJECT *p) {
      return p == Py_None;
    }

    GUESTOBJECT* full(const DM& m, bool simplify) {
      PyObject *p = swig::from(m);
      PyObject *method_name = PyUnicode_FromString("toarray");
      PyObject *cr = PyObject_CallMethodObjArgs(p, method_name, (simplify? Py_True: Py_False), 0);
      Py_DECREF(method_name);
      Py_DECREF(p);
      if (cr) return cr;
      return Py_None;
    }

    // Convert to a sparse matrix
    GUESTOBJECT* sparse(const DM& m) {
      PyObject *p = swig::from(m);
      PyObject *cr = PyObject_CallMethod(p, (char*) "tocsc", 0);
      Py_DECREF(p);
      if (cr) return cr;
      return Py_None;
    }

    PyObject* get_Python_helper(const std::string& name) {
%#if PY_VERSION_HEX < 0x03070000
      PyObject* module = PyImport_AddModule("casadi");
%#else
      PyObject* c_name = PyUnicode_FromString("casadi");
      PyObject* module = PyImport_GetModule(c_name);
      Py_DECREF(c_name);
%#endif
      if (!module) {
        if (PyErr_Occurred()) {
          PyErr_Clear();
        }
      }
      PyObject* dict = PyModule_GetDict(module);
      return PyDict_GetItemString(dict, (char*) name.c_str());
    }

    template<class T>
    bool casadi_object_from_fun(GUESTOBJECT *p, T** m, const std::string& fun, const std::function<bool(PyObject*, T**)> & conv) {
      PyObject* dm = get_Python_helper(fun);
      if (!dm) return false;
      PyObject *check_only = m? Py_False : Py_True;
      PyObject *cr = PyObject_CallFunctionObjArgs(dm, p, check_only, NULL);
      if (!cr) return false;
      bool ret;
      // None means not handled (#4216)
      if (cr == Py_None) {
        ret = false;
      } else if (PyBool_Check(cr)) {
        ret = PyObject_IsTrue(cr);
      } else {
        ret = conv(cr, m);
      }
      Py_DECREF(cr);
      return ret;
    }

    // False unless bytes holds exactly n values of type T
    template<typename T>
    bool from_bytes(PyObject* bytes, casadi_int n, std::vector<T>* data) {
      if (PyBytes_Size(bytes) != n * static_cast<casadi_int>(sizeof(T))) return false;
      data->resize(n);
      if (n) std::memcpy(data->data(), PyBytes_AsString(bytes), n * sizeof(T));
      return true;
    }

    bool SX_from_array_conv(GUESTOBJECT *p, casadi::SX** m) {
      std::vector<SXElem> data;
      if (!swig::assign(PyTuple_GetItem(p, 2), &data)) return false;
      casadi_int nrow; swig::assign(PyTuple_GetItem(p, 0), &nrow);
      casadi_int ncol; swig::assign(PyTuple_GetItem(p, 1), &ncol);
      if (m) {
        **m = casadi::SX::zeros(nrow, ncol);
        casadi_densify(get_ptr(data), (**m).sparsity().T(), (**m).ptr(), true);
      }
      return true;
    }

    bool IM_from_array_conv(GUESTOBJECT *p, casadi::IM** m) {
      if (!m) return true;
      casadi_int nrow; swig::assign(PyTuple_GetItem(p, 0), &nrow);
      casadi_int ncol; swig::assign(PyTuple_GetItem(p, 1), &ncol);
      std::vector<casadi_int> data;
      if (!from_bytes(PyTuple_GetItem(p, 3), nrow * ncol, &data)
          && !swig::assign(PyTuple_GetItem(p, 2), &data)) return false;
      **m = IM::zeros(nrow, ncol);
      casadi_densify(get_ptr(data), (**m).sparsity().T(), (**m).ptr(), true);
      return true;
    }

    bool DM_from_array_conv(GUESTOBJECT *p, casadi::DM** m) {
      if (!m) return true;
      casadi_int nrow; swig::assign(PyTuple_GetItem(p, 0), &nrow);
      casadi_int ncol; swig::assign(PyTuple_GetItem(p, 1), &ncol);
      std::vector<double> data;
      if (!from_bytes(PyTuple_GetItem(p, 2), nrow * ncol, &data)) return false;
      **m = DM::zeros(nrow, ncol);
      casadi_densify(get_ptr(data), (**m).sparsity().T(), (**m).ptr(), true);
      return true;
    }

    bool DM_from_csc_conv(GUESTOBJECT *p, casadi::DM** m) {
      std::vector<double> data;
      std::vector<casadi_int> colind, row;
      if (!swig::assign(PyTuple_GetItem(p, 4), &data)) return false;
      if (!swig::assign(PyTuple_GetItem(p, 3), &row)) return false;
      if (!swig::assign(PyTuple_GetItem(p, 2), &colind)) return false;
      casadi_int nrow; swig::assign(PyTuple_GetItem(p, 0), &nrow);
      casadi_int ncol; swig::assign(PyTuple_GetItem(p, 1), &ncol);
      **m = casadi::Matrix<double>(casadi::Sparsity(nrow,ncol,colind,row), data, false);
      return true;
    }

    bool SX_from_array(GUESTOBJECT *p, casadi::SX** m) {
      return casadi_object_from_fun<casadi::SX>(p, m, "SX_from_array", SX_from_array_conv);
    }

    bool IM_from_array(GUESTOBJECT *p, casadi::IM** m) {
      return casadi_object_from_fun<casadi::IM>(p, m, "IM_from_array", IM_from_array_conv);
    }

    bool DM_from_array(GUESTOBJECT *p, casadi::DM** m) {
      return casadi_object_from_fun<casadi::DM>(p, m, "DM_from_array", DM_from_array_conv);
    }

    bool DM_from_csc(GUESTOBJECT *p, casadi::DM** m) {
      return casadi_object_from_fun<casadi::DM>(p, m, "DM_from_csc", DM_from_csc_conv);
    }
  } // namespace casadi

  namespace swig {
    // A failed conversion leaves no Python error set
    inline bool casadi_ok(bool ok) {
      if (!ok && PyErr_Occurred()) PyErr_Clear();
      return ok;
    }

    // A float or int, but not a subclass, that fits in a double
    inline bool casadi_exact_scalar(GUESTOBJECT* p, double* d) {
      if (!PyFloat_CheckExact(p) && !PyLong_CheckExact(p)) return false;
      return casadi_ok(SWIG_IsOK(SWIG_AsVal(double)(p, d)));
    }

    // An exact float, int, DM, SX or MX converts or not by its type alone
    inline int casadi_host_conv_early(GUESTOBJECT* p, casadi::DM* m) {
      double d;
      if (casadi_exact_scalar(p, &d)) {
        if (m) *m = d;
        return true;
      }
      swig_type_info* t = casadi::exact_type(p);
      return t && t != $descriptor(casadi::Matrix<double>*) ? casadi_ok(false) : -1;
    }

    inline int casadi_host_conv_early(GUESTOBJECT* p, casadi::SX* m) {
      double d;
      if (casadi_exact_scalar(p, &d)) {
        if (m) *m = casadi::DM(d);
        return true;
      }
      return casadi::exact_type(p) == $descriptor(casadi::MX*) ? casadi_ok(false) : -1;
    }

    inline int casadi_host_conv_early(GUESTOBJECT* p, casadi::MX* m) {
      double d;
      if (casadi_exact_scalar(p, &d)) {
        if (m) *m = casadi::DM(d);
        return true;
      }
      return casadi::exact_type(p) == $descriptor(casadi::Matrix<casadi::SXElem>*) ? casadi_ok(false) : -1;
    }

    inline int casadi_host_conv_early(GUESTOBJECT* p, casadi::Sparsity* m) {
      return casadi::exact_type(p) ? casadi_ok(false) : -1;
    }

    inline int casadi_host_conv(GUESTOBJECT* p, casadi::Slice* m) {
      // Python casadi_int
      if (PyLong_Check(p)) {
        casadi_int i = PyLong_AsLongLong(p);
        if (i == -1 && PyErr_Occurred()) return casadi_ok(false);
        if (m) {
          m->start = i;
          m->stop = m->start+1;
          if (m->stop==0) m->stop = std::numeric_limits<casadi_int>::max();
        }
        return true;
      }
      // Python slice - use Limited API compatible approach
      if (PySlice_Check(p)) {
        Py_ssize_t start, stop, step;
%#if PY_VERSION_HEX >= 0x03060100
        int res = PySlice_Unpack(p, &start, &stop, &step);
%#else
        // Python 2.7 and early Python 3.x use _PySlice_Unpack (private API)
        int res = _PySlice_Unpack(p, &start, &stop, &step);
%#endif
        if (res < 0) {
          return false;  // TypeError already set by PySlice_Unpack
        }

        if (m) {
          // Map PySlice_Unpack sentinels (either extreme) to CasADi's limits
          m->start = (start == PY_SSIZE_T_MIN || start == PY_SSIZE_T_MAX)
              ? std::numeric_limits<casadi_int>::min()
              : static_cast<casadi_int>(start);
          m->stop = (stop == PY_SSIZE_T_MAX || stop == PY_SSIZE_T_MIN)
              ? std::numeric_limits<casadi_int>::max()
              : static_cast<casadi_int>(stop);
          if (step != 1) {
            m->step = static_cast<casadi_int>(step);
          }
        }
        return true;
      }
      return false;
    }

    inline int casadi_host_conv(GUESTOBJECT* p, casadi::DM* m) {
      using namespace casadi;
      // Object has __DM__ method; None means not a DM
      if (PyErr_Occurred()) PyErr_Clear(); // Clear pending exception before type check
      if (PyObject_HasAttrString(p,"__DM__")) {
        char name[] = "__DM__";
        PyObject *cr = PyObject_CallMethod(p, name, 0);
        if (!cr) return false;
        if (cr == Py_None) { Py_DECREF(cr); return false; }
        bool result = swig::assign(cr, m);
        Py_DECREF(cr);
        return result;
      }

      DM* mp = m;
      if (DM_from_array(p, m ? &mp : 0)) return true;

      if (DM_from_csc(p, m ? &mp : 0)) return true;

      {
        std::vector<double> t;
        bool res = swig::assign(p, m ? &t : 0);
        if (m) *m = t.empty() ? casadi::Matrix<double>(0,0) : casadi::Matrix<double>(t);
        return res;
      }
    }

    inline int casadi_host_conv(GUESTOBJECT* p, casadi::IM* m) {
      using namespace casadi;
      // Numpy arrays will be cast to dense Matrix<casadi_int>
      IM* mp = m;
      if (IM_from_array(p, m ? &mp : 0)) return true;

      if (PyErr_Occurred()) PyErr_Clear(); // Clear pending exception before type check
      if (PyObject_HasAttrString(p,"__IM__")) {
        PyObject *cr = PyObject_CallMethod(p, (char*) "__IM__", 0);
        if (!cr) return false;
        bool result = swig::assign(cr, m);
        Py_DECREF(cr);
        return result;
      }

      {
        std::vector <casadi_int> t;
        if (swig::assign(p, m ? &t : 0)) {
          if (m) *m = casadi::Matrix<casadi_int>(t);
          return true;
        }
      }

      return -1;
    }

    inline int casadi_host_conv(GUESTOBJECT* p, casadi::SX* m) {
      using namespace casadi;
      // Object has __SX__ method; checked before the numeric heuristics
      if (PyErr_Occurred()) PyErr_Clear(); // Clear pending exception before type check
      if (PyObject_HasAttrString(p,"__SX__")) {
        PyObject *cr = PyObject_CallMethod(p, (char*) "__SX__", 0);
        if (!cr) return false;
        if (cr == Py_None) { Py_DECREF(cr); }
        else {
          // swig::assign COPIES into *m; the hook result dies with cr
          bool flag = swig::assign(cr, m);
          Py_DECREF(cr);
          if (flag) return true;
        }
      }

      return -1;
    }

    inline int casadi_host_conv_late(GUESTOBJECT* p, casadi::SX* m) {
      using namespace casadi;
      // Clear any error from a failed DM conversion
      if (PyErr_Occurred()) PyErr_Clear();
      // Numpy arrays will be cast to dense SX
      SX* mp = m;
      if (SX_from_array(p, m ? &mp : 0)) return true;

      return -1;
    }

    inline int casadi_host_conv(GUESTOBJECT* p, casadi::MX* m) {
      using namespace casadi;
      // Object has __MX__ method; checked before the numeric heuristics
      if (PyErr_Occurred()) PyErr_Clear(); // Clear pending exception before type check
      if (PyObject_HasAttrString(p,"__MX__")) {
        PyObject *cr = PyObject_CallMethod(p, (char*) "__MX__", 0);
        if (!cr) return false;
        if (cr == Py_None) { Py_DECREF(cr); }
        else {
          // swig::assign copies into *m
          bool flag = swig::assign(cr, m);
          Py_DECREF(cr);
          if (flag) return true;
        }
      }

      return -1;
    }

    inline int casadi_host_conv_late(GUESTOBJECT* p, casadi::MX* m) {
      // Clear any error from a failed DM conversion
      if (PyErr_Occurred()) PyErr_Clear();

      return -1;
    }

    inline GUESTOBJECT* casadi_host_from(const casadi::GenericType& a) {
      using namespace casadi;
      switch (a.getType()) {
      case OT_NULL:
      case OT_VOIDPTR:
        // A new reference, like every other branch
        Py_INCREF(Py_None);
        return Py_None;
      default: return 0;
      }
    }
  } // namespace swig
}

// Indexing helpers of the matrix classes
%define %matrix_helpers(Type, xIndex)
%pythoncode %{
    @property
    def shape(self):
        return (self.size1(),self.size2())

    def reshape(self,arg):
        return _casadi.reshape(self,arg)

    @property
    def T(self):
        return _casadi.transpose(self)

    def __getitem__(self, s):
          if isinstance(s, tuple) and len(s)==2:
            if s[1] is None: raise TypeError("Cannot slice with None")
            return self.get(False, s[0], s[1])
          return self.get(False, s)

    def __iter__(self):
      raise TypeError("""CasADi matrices are not iterable by design.
                      Did you mean to iterate over m.nz, with m IM/DM/SX?
                      Did you mean to iterate over horzsplit(m,1)/vertsplit(m,1) with m IM/DM/SX/MX?
                      Did you expect numpy behaviour (looping over rows)? Use ca.array().
                      """)

    def __setitem__(self,s,val):
          if isinstance(s,tuple) and len(s)==2:
            return self.set(val, False, s[0], s[1])
          return self.set(val, False, s)

    @property
    def nz(self):
      return NZproxy(self)

%}
%stub_attr(shape, tuple[int, int])
%stub_attr(T,     Self)
%stub_attr(nz,    NZproxy[Self])
%stub_method(reshape,     Self, arg: tuple[int, int] | int | Sparsity)
// Indexing; MX adds MX indices in its own %extend
%stub_method(__getitem__, Self, s: xIndex | tuple[xIndex, xIndex])
%stub_method(__setitem__, None, s: xIndex | tuple[xIndex, xIndex], val: _DM | SX | MX)
%stub_method0(__iter__, %arg(Iterator[Self]))
%enddef

// Zero-argument getter as a property
%define make_property(class, name, type)
  %rename(_ ## name) class ## :: ## name;
  %extend class {
    %pythoncode %{
      @property
      def name(self):
        return self._ ## name()
    %}
    %stub_attr(name, type)
  }
%enddef
