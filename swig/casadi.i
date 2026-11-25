/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


%module(package="casadi",directors=1) casadi

// Per-target parts: head (here), pre (before stl.i), decl (before the first class), post (after the last)
#ifdef SWIGPYTHON
%include "python/casadi_python_head.i"
#elif defined(SWIGMATLAB)
%include "matlab/casadi_matlab_head.i"
#elif defined(SWIGJULIA)
%include "julia/casadi_julia_head.i"
#elif defined(SWIGWASMJS)
%include "wasm-js/casadi_wasmjs_head.i"
#elif defined(SWIGJAVA)
%include "java/casadi_java_head.i"
#elif defined(SWIGRUST)
%include "rust/casadi_rust_head.i"
#elif defined(SWIGDOTNET)
%include "dotnet/casadi_dotnet_head.i"
#endif


#ifdef CASADI_WITH_COPYSIGN_UNDEF
%{
#ifdef copysign
#undef copysign
#endif
%}
#endif // CASADI_WITH_COPYSIGN_UNDEF

// Include all public CasADi C++
%{
#include <casadi/casadi.hpp>
#include <casadi/core/casadi_interrupt.hpp>
%}

// casadi_int type
%include <casadi/core/casadi_types.hpp>

// No-op %stub_* macros without -stubs
#ifndef SWIG_STUBS_ENABLED
%define %stub_method(NAME, RET, ARGS...) %enddef
%define %stub_method0(NAME, RET) %enddef
%define %stub_overload_method(NAME, RET, ARGS...) %enddef
%define %stub_overload_method0(NAME, RET) %enddef
%define %stub_overload_method_selftyped(NAME, RET, SELFTYPE, ARGS...) %enddef
%define %stub_func(NAME, RET, ARGS...) %enddef
%define %stub_overload_func(NAME, RET, ARGS...) %enddef
%define %stub_attr(NAME, TYPE...) %enddef
%define %stub_const(NAME, TYPE...) %enddef
%define %stub_class_begin(NAME) %enddef
%define %stub_alias_in(NAME, TYPES...) %enddef
%define %stub_alias_out(NAME, TYPES...) %enddef
%define %stub_alias_in_key(KEY, TYPES...) %enddef
%define %stub_alias_out_key(KEY, TYPES...) %enddef
%define %ts_alias_in(NAME) %enddef
%define %ts_alias_out(NAME) %enddef
#endif

  /// Data structure in the target language holding data; targets define it in their head file
#ifndef GUESTOBJECT
#define GUESTOBJECT void
#endif
#ifndef CASADI_OK
#define CASADI_OK(x) x
#endif

// Define printing routine


// Spurious warnings, e.g. for vertcat(SXVector)
#pragma SWIG nowarn=509,303,302

#define CASADI_EXPORT

// Incude cmath early on, see #622
%begin %{
#include <cmath>
#ifdef _XOPEN_SOURCE
#undef _XOPEN_SOURCE
#endif
#ifdef _POSIX_C_SOURCE
#undef _POSIX_C_SOURCE
#endif
%}

%ignore *::operator->;

%rename(str) get_str;


// Doc styles: error, overview, group

%feature("autodoc", "1");

%feature("customdoc", "1");

// Docstrings keep the \lang groups of the target language; pre files override
%feature("docstring:lang", "python");


%feature("customdoc:arg:self", "self");

#ifdef CASADI_MNAME
  #define MNAME CASADI_MNAME
#else
  #define MNAME "$name"
  %feature("customdoc:proto:constructor", "$name($in)");
  %feature("customdoc:proto:single_out", "$name($in) -> $out");
  %feature("customdoc:proto:normal", "$name($in) -> ($out)");
  %feature("customdoc:main", "  $brief\n\n::\n\n$overview\n$main");
#endif

%feature("customdoc:arg:normal:style_error", "$type");
%feature("customdoc:arg:only:style_error", "$type");
%feature("customdoc:arg:separator:style_error", ",");
%feature("customdoc:proto:void:style_error", MNAME "($in)");
%feature("customdoc:proto:single_out:style_error", MNAME "($in)");
%feature("customdoc:proto:normal:style_error", MNAME "($in)");
%feature("customdoc:proto:constructor:style_error", MNAME "($in)");

%feature("customdoc:arg:normal", "$type $name");
%feature("customdoc:arg:only", "$type $name");
%feature("customdoc:arg:only:out", "$type");
%feature("customdoc:arg:no_name", "out$ip");
%feature("customdoc:arg:separator", ", ");


%feature("customdoc:proto:void", MNAME "($in)");
%feature("customdoc:proto:single_out:style_group", MNAME "($in)");
%feature("customdoc:proto:normal:style_group", MNAME "($in)");
%feature("customdoc:proto:constructor:style_group", MNAME "($in)");

%feature("customdoc:protoline", "    $proto");
%feature("customdoc:protoline:style_overview", "  $proto");
%feature("customdoc:protoline:nobrief:style_overview", "  $proto");

%feature("customdoc:protoline:style_group", "  $proto");

%feature("customdoc:group", "\n.......\n\n::\n\n$group\n$main\n\n.............\n\n");

// append works for strings

%naturalvar;

// Make data members read-only
%immutable;

// Make sure that a copy constructor is created
%copyctor;

#ifndef SWIGXML
%feature("compactdefaultargs","1");
#endif //SWIGXML

// STL
#ifdef SWIGXML
namespace std {
  template<class T> class vector {};
  template<class A, class B> class pair {};
  template<class A, class B> class map {};
}
#else // SWIGXML
#ifdef SWIGPYTHON
%include "python/casadi_python_pre.i"
#elif defined(SWIGMATLAB)
%include "matlab/casadi_matlab_pre.i"
#elif defined(SWIGJULIA)
%include "julia/casadi_julia_pre.i"
#elif defined(SWIGWASMJS)
%include "wasm-js/casadi_wasmjs_pre.i"
#elif defined(SWIGJAVA)
%include "java/casadi_java_pre.i"
#elif defined(SWIGRUST)
%include "rust/casadi_rust_pre.i"
#elif defined(SWIGDOTNET)
%include "dotnet/casadi_dotnet_pre.i"
#endif
%include "stl.i"
#endif // SWIGXML

%include "doc.i"


// Exceptions handling
%include "exception.i"
%exception {
  try {
    $action
  } catch(const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

// Python expects __int__ to raise on failure
%exception __int__ {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

// numpy printing only catches TypeError/ValueError (#701)
%exception __bool__ {
 try {
    $action
  } catch (const std::exception& e) {
   SWIG_exception(SWIG_TypeError, e.what());
  }
}


#ifndef SWIGXML

// Can be overloaded by specifying before importing casadi.i
%fragment("casadi_extra_decl", "header") {}

%fragment("casadi_decl", "header", fragment="casadi_extra_decl", fragment="StdTraits") {
  namespace swig {
    // *m is the caller's object filled in place, or redirected to the wrapped object
    template<class T, class C = typename traits<T>::category, class = void> struct traits_asref {
      static bool asref(GUESTOBJECT *p, T** m) {
        T* r = 0;
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
    template<class T> struct traits_asref<T, value_category, void> {
      static bool asref(GUESTOBJECT *p, T** m) { return SWIG_IsOK(asval(p, *m)); }
    };
    template<class T> struct traits_asref<T, pointer_category, decltype((void)&traits_asptr<T>::in_place)> {
      static bool asref(GUESTOBJECT *p, T** m) { return traits_asptr<T>::in_place(p, m); }
    };

    // Convert p via *m (see traits_asref); m==0: check only
    template<class T> inline bool asref(GUESTOBJECT *p, T** m) {
      return m ? traits_asref<T>::asref(p, m) : check<T>(p);
    }

    // Convert p into *v; v==0: check only
    template<class T> inline bool assign(GUESTOBJECT *p, T* v) {
      T* q = v;
      if (!asref(p, v ? &q : 0)) return false;
      if (q != v) *v = *q;
      return true;
    }
  } // namespace swig

  namespace casadi {
    /* Check if Null or None */
    bool is_null(GUESTOBJECT *p);

    GUESTOBJECT* full(const DM& m, bool simplify=false);

    // Convert to a sparse matrix
    GUESTOBJECT* sparse(const DM& m);

    inline GUESTOBJECT* full_or_sparse(const DM& m, bool simplify=false) {
      if (m.is_dense()) {
        return full(m, simplify);
      } else {
        GUESTOBJECT* p = sparse(m);
        if (is_null(p)) return swig::from(m);
        return p;
      }
    }
  } // namespace casadi
 }

// Conversions: swig::traits specializations for casadi's typemaps

%fragment("casadi_aux", "header", fragment="casadi_host", fragment="StdTraits") {
  namespace swig {
    // Proxy pointer (SWIG_OLDOBJ) or a new object filled by conv (SWIG_NEWOBJ); m==0: check only
    template<class M> int casadi_asptr(GUESTOBJECT* p, M** m, swig_type_info* desc, bool (*conv)(GUESTOBJECT*, M*)) {
      if (casadi::is_null(p)) return SWIG_ERROR;
      if (desc) {
        void* vp = 0;
        if (SWIG_IsOK(SWIG_ConvertPtr(p, m ? &vp : 0, desc, 0))) {
          if (m) *m = static_cast<M*>(vp);
          return SWIG_OLDOBJ;
        }
      }
      if (!conv) return SWIG_ERROR;
      if (!m) return CASADI_OK(conv(p, 0)) ? SWIG_OK : SWIG_ERROR;
      M* r = new M();
      if (!CASADI_OK(conv(p, r))) {
        delete r;
        return SWIG_ERROR;
      }
      *m = r;
      return SWIG_NEWOBJ;
    }

    // Same, with asref semantics: *m is redirected to the proxy or conv fills **m in place
    template<class M> bool casadi_in_place(GUESTOBJECT* p, M** m, swig_type_info* desc, bool (*conv)(GUESTOBJECT*, M*)) {
      if (casadi::is_null(p)) return false;
      if (desc) {
        void* vp = 0;
        if (SWIG_IsOK(SWIG_ConvertPtr(p, &vp, desc, 0))) {
          *m = static_cast<M*>(vp);
          return true;
        }
      }
      return conv && CASADI_OK(conv(p, *m));
    }

    // Host conversions, overloaded per target in casadi_host: 1 converted, 0 rejected, -1 not handled
    template<class M> int casadi_host_conv(GUESTOBJECT*, M*) { return -1; }
    template<class M> int casadi_host_conv_late(GUESTOBJECT*, M*) { return -1; }
    template<class M> GUESTOBJECT* casadi_host_from(const M&) { return 0; }
  } // namespace swig
}

// Container policy: casadi's acceptance rules on top of SWIG's generic sequence/map/pair traits
%fragment("casadi_map_traits", "header", fragment="StdMapTraits") {}
%fragment("casadi_containers", "header", fragment="casadi_aux", fragment="StdVectorTraits",
          fragment="casadi_map_traits", fragment="StdPairTraits", fragment=SWIG_Traits_frag(bool),
          fragment=SWIG_Traits_frag(long long), fragment=SWIG_Traits_frag(double),
          fragment=SWIG_Traits_frag(std::string)) {
}

// Element traits, one fragment per type, named after the spelling in the %template() lines

// Wrapped classes; host conversions only (Slice)
%define %casadi_class_traits(xName, xType...)
%fragment(SWIG_Traits_frag(xType), "header", fragment="casadi_aux") {
  namespace swig {
    template<> struct traits<xType> {
      typedef pointer_category category;
      static const char* type_name() { return xName; }
    };
    template<> struct traits_asptr<xType> {
      static bool conv(GUESTOBJECT* p, xType* m) {
        return casadi_host_conv(p, m) > 0;
      }
      static int asptr(GUESTOBJECT* p, xType** m) {
        return casadi_asptr<xType>(p, m, $descriptor(xType*), conv);
      }
      static bool in_place(GUESTOBJECT* p, xType** m) {
        return casadi_in_place<xType>(p, m, $descriptor(xType*), conv);
      }
    };
    template<> struct traits_from<xType> {
      static GUESTOBJECT* from(const xType& a) {
        return SWIG_NewPointerObj(new xType(a), $descriptor(xType*), SWIG_POINTER_OWN);
      }
    };
  } // namespace swig
}
%enddef

%casadi_class_traits("Sparsity", casadi::Sparsity)
%casadi_class_traits("Function", casadi::Function)
%casadi_class_traits("Linsol", casadi::Linsol)
%casadi_class_traits("Slice", casadi::Slice)

%fragment(SWIG_Traits_frag(casadi::Matrix<double>), "header", fragment="casadi_aux",
          fragment=SWIG_Traits_frag(double), fragment="casadi_containers") {
  namespace swig {
    template<> struct traits<casadi::Matrix<double> > {
      typedef pointer_category category;
      static const char* type_name() { return "DM"; }
    };
    template<> struct traits_asptr<casadi::Matrix<double> > {
      static bool conv(GUESTOBJECT* p, casadi::DM* m) {
        using namespace casadi;
        // Object is a sparsity pattern
        {
          void *sp = 0;
          if (SWIG_IsOK(SWIG_ConvertPtr(p, &sp, $descriptor(casadi::Sparsity*), 0))) {
            if (m) *m = DM::ones(*static_cast<Sparsity*>(sp));
            return true;
          }
        }

        // Double scalar
        {
          double tmp;
          if (swig::assign(p, m? &tmp: 0)) {
            if (m) *m = tmp;
            return true;
          }
        }

        // Host representations
        return casadi_host_conv(p, m) > 0;
      }
      static int asptr(GUESTOBJECT* p, casadi::DM** m) {
        return casadi_asptr<casadi::DM>(p, m, $descriptor(casadi::Matrix<double>*), conv);
      }
      static bool in_place(GUESTOBJECT* p, casadi::DM** m) {
        return casadi_in_place<casadi::DM>(p, m, $descriptor(casadi::Matrix<double>*), conv);
      }
    };
    template<> struct traits_from<casadi::Matrix<double> > {
      static GUESTOBJECT* from(const casadi::DM& a) {
        return SWIG_NewPointerObj(new casadi::DM(a), $descriptor(casadi::Matrix<double>*), SWIG_POINTER_OWN);
      }
    };
  } // namespace swig
}

%fragment(SWIG_Traits_frag(casadi::Matrix<casadi_int>), "header", fragment="casadi_aux",
          fragment=SWIG_Traits_frag(long long), fragment=SWIG_Traits_frag(casadi::Matrix<double>)) {
  namespace swig {
    template<> struct traits<casadi::Matrix<casadi_int> > {
      typedef pointer_category category;
      static const char* type_name() { return "IM"; }
    };
    template<> struct traits_asptr<casadi::Matrix<casadi_int> > {
      static bool conv(GUESTOBJECT* p, casadi::IM* m) {
        using namespace casadi;
        // Object is a sparsity pattern
        {
          void *sp = 0;
          if (SWIG_IsOK(SWIG_ConvertPtr(p, &sp, $descriptor(casadi::Sparsity*), 0))) {
            if (m) *m = IM::ones(*static_cast<Sparsity*>(sp));
            return true;
          }
        }

        // First convert to integer
        {
          casadi_int tmp;
          if (swig::assign(p, m? &tmp: 0)) {
            if (m) *m = tmp;
            return true;
          }
        }

        // Host representations
        {
          int r = casadi_host_conv(p, m);
          if (r >= 0) return r;
        }

        // Convert from DM
        {
          DM tmp;
          if (swig::assign(p, m? &tmp: 0)) {
            // Check integrality
            for (double d : tmp.nonzeros()) {
              if (d!=casadi_int(d)) return false;
            }
            // Convert
            if (m) *m = casadi::Matrix<double>(tmp);
            return true;
          }
        }

        // No match
        return false;
      }
      static int asptr(GUESTOBJECT* p, casadi::IM** m) {
        return casadi_asptr<casadi::IM>(p, m, 0, conv);
      }
      static bool in_place(GUESTOBJECT* p, casadi::IM** m) {
        return casadi_in_place<casadi::IM>(p, m, 0, conv);
      }
    };
    // An IM leaves as a DM
    template<> struct traits_from<casadi::Matrix<casadi_int> > {
      static GUESTOBJECT* from(const casadi::IM& a) {
        return swig::from(casadi::DM(a));
      }
    };
  } // namespace swig
}

%fragment(SWIG_Traits_frag(casadi::Matrix<casadi::SXElem>), "header", fragment="casadi_aux",
          fragment=SWIG_Traits_frag(casadi::Matrix<double>)) {
  namespace swig {
    template<> struct traits<casadi::Matrix<casadi::SXElem> > {
      typedef pointer_category category;
      static const char* type_name() { return "SX"; }
    };
    template<> struct traits_asptr<casadi::Matrix<casadi::SXElem> > {
      static bool conv(GUESTOBJECT* p, casadi::SX* m) {
        // Host representations
        {
          int r = casadi_host_conv(p, m);
          if (r >= 0) return r;
        }

        // Try first converting to a temporary DM
        {
          casadi::DM tmp;
          if (swig::assign(p, m? &tmp: 0)) {
            if (m) *m = tmp;
            return true;
          }
        }
        {
          int r = casadi_host_conv_late(p, m);
          if (r >= 0) return r;
        }
        // No match
        return false;
      }
      static int asptr(GUESTOBJECT* p, casadi::SX** m) {
        return casadi_asptr<casadi::SX>(p, m, $descriptor(casadi::Matrix<casadi::SXElem>*), conv);
      }
      static bool in_place(GUESTOBJECT* p, casadi::SX** m) {
        return casadi_in_place<casadi::SX>(p, m, $descriptor(casadi::Matrix<casadi::SXElem>*), conv);
      }
    };
    template<> struct traits_from<casadi::Matrix<casadi::SXElem> > {
      static GUESTOBJECT* from(const casadi::SX& a) {
        return SWIG_NewPointerObj(new casadi::SX(a), $descriptor(casadi::Matrix<casadi::SXElem> *), SWIG_POINTER_OWN);
      }
    };
  } // namespace swig
}

%fragment(SWIG_Traits_frag(casadi::SXElem), "header", fragment="casadi_aux",
          fragment=SWIG_Traits_frag(casadi::Matrix<casadi::SXElem>)) {
  namespace swig {
    template<> struct traits<casadi::SXElem> {
      typedef pointer_category category;
      static const char* type_name() { return "SXElem"; }
    };
    template<> struct traits_asptr<casadi::SXElem> {
      static bool conv(GUESTOBJECT* p, casadi::SXElem* m) {
        // Try first converting to a temporary SX
        casadi::SX tmp, *mt=&tmp;
        if(swig::asref(p, m ? &mt : 0)) {
          if (m && !mt->is_scalar()) return false;
          if (m) *m = mt->scalar();
          return true;
        }
        // No match
        return false;
      }
      static int asptr(GUESTOBJECT* p, casadi::SXElem** m) {
        return casadi_asptr<casadi::SXElem>(p, m, 0, conv);
      }
      static bool in_place(GUESTOBJECT* p, casadi::SXElem** m) {
        return casadi_in_place<casadi::SXElem>(p, m, 0, conv);
      }
    };
    // An SXElem leaves as an SX
    template<> struct traits_from<casadi::SXElem> {
      static GUESTOBJECT* from(const casadi::SXElem& a) {
        return swig::from(casadi::SX(a));
      }
    };
  } // namespace swig
}

%fragment(SWIG_Traits_frag(casadi::MX), "header", fragment="casadi_aux",
          fragment=SWIG_Traits_frag(casadi::Matrix<double>)) {
  namespace swig {
    template<> struct traits<casadi::MX> {
      typedef pointer_category category;
      static const char* type_name() { return "MX"; }
    };
    template<> struct traits_asptr<casadi::MX> {
      static bool conv(GUESTOBJECT* p, casadi::MX* m) {
        // Host representations
        {
          int r = casadi_host_conv(p, m);
          if (r >= 0) return r;
        }

        // Try first converting to a temporary DM
        {
          casadi::DM tmp;
          if (swig::assign(p, m ? &tmp : 0)) {
            if (m) *m = tmp;
            return true;
          }
        }
        {
          int r = casadi_host_conv_late(p, m);
          if (r >= 0) return r;
        }
        // No match
        return false;
      }
      static int asptr(GUESTOBJECT* p, casadi::MX** m) {
        return casadi_asptr<casadi::MX>(p, m, $descriptor(casadi::MX*), conv);
      }
      static bool in_place(GUESTOBJECT* p, casadi::MX** m) {
        return casadi_in_place<casadi::MX>(p, m, $descriptor(casadi::MX*), conv);
      }
    };
    template<> struct traits_from<casadi::MX> {
      static GUESTOBJECT* from(const casadi::MX& a) {
        return SWIG_NewPointerObj(new casadi::MX(a), $descriptor(casadi::MX*), SWIG_POINTER_OWN);
      }
    };
  } // namespace swig
}

// GenericType: declarations here, bodies in casadi_generictype (they need the container traits of GenericType itself)
%fragment(SWIG_Traits_frag(casadi::GenericType), "header", fragment="casadi_aux") {
  namespace swig {
    template<> struct traits<casadi::GenericType> {
      typedef pointer_category category;
      static const char* type_name() { return "GenericType"; }
    };
    template<> struct traits_asptr<casadi::GenericType> {
      static bool conv(GUESTOBJECT* p, casadi::GenericType* m);
      static int asptr(GUESTOBJECT* p, casadi::GenericType** m) {
#ifdef CASADI_NOT_NONE
        if (CASADI_NOT_NONE(p)) return casadi_asptr<casadi::GenericType>(p, m, $descriptor(casadi::GenericType*), conv);
        if (m) *m = new casadi::GenericType();
        return m ? SWIG_NEWOBJ : SWIG_OK;
#else
        return casadi_asptr<casadi::GenericType>(p, m, $descriptor(casadi::GenericType*), conv);
#endif
      }
      static bool in_place(GUESTOBJECT* p, casadi::GenericType** m) {
#ifdef CASADI_NOT_NONE
        if (CASADI_NOT_NONE(p)) return casadi_in_place<casadi::GenericType>(p, m, $descriptor(casadi::GenericType*), conv);
        **m = casadi::GenericType();
        return true;
#else
        return casadi_in_place<casadi::GenericType>(p, m, $descriptor(casadi::GenericType*), conv);
#endif
      }
    };
    template<> struct traits_from<casadi::GenericType> {
      static GUESTOBJECT* from(const casadi::GenericType& a);
    };
  } // namespace swig
}

%fragment("casadi_generictype", "header", fragment=SWIG_Traits_frag(casadi::GenericType),
          fragment="casadi_containers", fragment=SWIG_Traits_frag(bool), fragment=SWIG_Traits_frag(long long),
          fragment=SWIG_Traits_frag(double), fragment=SWIG_Traits_frag(std::string),
          fragment=SWIG_Traits_frag(casadi::Function)) {
  namespace swig {
    // Try one alternative of GenericType
    template<class T> bool casadi_generic(GUESTOBJECT* p, casadi::GenericType* m) {
#ifdef SWIGMATLAB
      if (!m) return swig::check<T>(p);
      T tmp;
      if (!SWIG_IsOK(swig::asval<T>(p, &tmp))) return false;
      *m = casadi::GenericType(tmp);
#else
      if (!m) return CASADI_OK(swig::check<T>(p));
      T tmp, *tmp_ptr = &tmp;
      if (!CASADI_OK(swig::asref(p, &tmp_ptr))) return false;
      *m = casadi::GenericType(*tmp_ptr);
#endif
      return true;
    }

    bool traits_asptr<casadi::GenericType>::conv(GUESTOBJECT* p, casadi::GenericType* m) {
      using casadi::GenericType;
      {
        int r = casadi_host_conv(p, m);
        if (r >= 0) return r;
      }
      // Try to convert to different types
      return casadi_generic<bool>(p, m)
          || casadi_generic<casadi_int>(p, m)
          || casadi_generic<double>(p, m)
          || casadi_generic<std::string>(p, m)
#ifdef SWIGPYTHON
          || casadi_generic<GenericType::Dict>(p, m) // a dict before any sequence
#endif
          || casadi_generic<std::vector<bool> >(p, m)
          || casadi_generic<std::vector<casadi_int> >(p, m)
          || casadi_generic<std::vector<double> >(p, m)
          || casadi_generic<std::vector<std::string> >(p, m)
          || casadi_generic<std::vector<std::vector<std::string> > >(p, m)
          || casadi_generic<std::vector<std::vector<casadi_int> > >(p, m)
          || casadi_generic<std::vector<std::vector<double> > >(p, m)
          || casadi_generic<casadi::Function>(p, m)
          || casadi_generic<std::vector<casadi::Function> >(p, m)
#ifndef SWIGPYTHON
          || casadi_generic<GenericType::Dict>(p, m)
#endif
          || casadi_generic<std::vector<GenericType::Dict> >(p, m)
          || casadi_generic<std::vector<std::vector<GenericType> > >(p, m)
          || casadi_generic<std::vector<GenericType> >(p, m);
    }

    GUESTOBJECT* traits_from<casadi::GenericType>::from(const casadi::GenericType& a) {
      using namespace casadi;
      switch (a.getType()) {
      case OT_BOOL: return swig::from(a.as_bool());
      case OT_INT: return swig::from(a.as_int());
      case OT_DOUBLE: return swig::from(a.as_double());
      case OT_STRING: return swig::from(a.as_string());
      case OT_INTVECTOR: return swig::from(a.as_int_vector());
      case OT_INTVECTORVECTOR: return swig::from(a.as_int_vector_vector());
      case OT_BOOLVECTOR: return swig::from(a.as_bool_vector());
      case OT_DOUBLEVECTOR: return swig::from(a.as_double_vector());
      case OT_DOUBLEVECTORVECTOR: return swig::from(a.as_double_vector_vector());
      case OT_STRINGVECTOR: return swig::from(a.as_string_vector());
      case OT_STRINGVECTORVECTOR: return swig::from(a.as_string_vector_vector());
      case OT_DICT: return swig::from(a.as_dict());
      case OT_DICTVECTOR: return swig::from(a.as_dict_vector());
      case OT_VECTORVECTOR: return swig::from(a.as_vector_vector());
      case OT_VECTOR: return swig::from(a.as_vector());
      case OT_FUNCTION: return swig::from(a.as_function());
      case OT_FUNCTIONVECTOR: return swig::from(a.as_function_vector());
      default: return casadi_host_from(a);
      }
    }
  } // namespace swig
}

// Collect all fragments
%fragment("casadi_all", "header", fragment="casadi_aux,casadi_extra,casadi_containers,casadi_generictype",
          fragment=SWIG_Traits_frag(casadi::Sparsity), fragment=SWIG_Traits_frag(casadi::Function),
          fragment=SWIG_Traits_frag(casadi::Slice), fragment=SWIG_Traits_frag(casadi::Matrix<double>),
          fragment=SWIG_Traits_frag(casadi::Matrix<casadi_int>), fragment=SWIG_Traits_frag(casadi::Matrix<casadi::SXElem>),
          fragment=SWIG_Traits_frag(casadi::SXElem), fragment=SWIG_Traits_frag(casadi::MX)
          ) { }

// Can be overloaded by specifying before importing casadi.i
%fragment("casadi_extra", "header") {}


#endif // SWIGXML

 // Input typemaps; xStubIn/xTsStub are the .pyi/.d.ts annotations
%define %casadi_input_typemaps(xName, xStubIn, xTsStub, xPrec, xType...)
#ifdef CASADI_TARGET_INPUT_TYPEMAPS
CASADI_TARGET_INPUT_TYPEMAPS(xName, xStubIn, xTsStub, xType)
#endif
 // Pass input by value, check if matches
%typemap(typecheck, noblock=1, precedence=xPrec, fragment="casadi_all") xType {
  $1 = swig::check< xType >($input);
 }

 // Directorout typemap; as input by value
%typemap(directorout, noblock=1, fragment="casadi_all") xType {
    if (!swig::assign< xType >($input, &$result)) {
      %dirout_fail(SWIG_TypeError,"$type");
    }
 }

 // Pass input by value, convert argument
%typemap(in, doc=xName, pystub_in=xStubIn, tsstub_in=xTsStub, noblock=1, fragment="casadi_all") xType {
  if (!swig::assign< xType >($input, &$1)) SWIG_exception_fail(SWIG_TypeError,"Failed to convert input $argnum to type '" xName "'.");
 }

 // Pass input by value, cleanup
%typemap(freearg, noblock=1) xType {}

 // Pass input by reference, check if matches
%typemap(typecheck, noblock=1, precedence=xPrec, fragment="casadi_all") const xType& {
  $1 = swig::check< xType >($input);
 }

 // Pass input by reference, convert argument
%typemap(in, doc=xName, pystub_in=xStubIn, tsstub_in=xTsStub, noblock=1, fragment="casadi_all") const xType & (xType m) {
  $1 = &m;
  if (!swig::asref< xType >($input, &$1)) SWIG_exception_fail(SWIG_TypeError,"Failed to convert input $argnum to type '" xName "'.");
 }

 // Pass input by reference, cleanup
%typemap(freearg, noblock=1) const xType & {}

%enddef

 // Output typemaps; xStubOut/xTsStub are the .pyi/.d.ts annotations
%define %casadi_output_typemaps(xName, xStubIn, xStubOut, xTsStub, xType...)


// Inputs marked OUTPUT are also returned by the function, ...
%typemap(argout, noblock=1,fragment="casadi_all") xType &OUTPUT {
  %append_output(swig::from(*$1));
 }

// ... and the corresponding inputs are ignored.
%typemap(in, doc=xName, pystub_out=xStubOut, tsstub_out=xTsStub, noblock=1, numinputs=0) xType &OUTPUT (xType m) {
 $1 = &m;
}

 // Directorin typemap; as output
%typemap(directorin, noblock=1, fragment="casadi_all") xType, const xType {
    if(!($input = swig::from< xType >($1))) %dirout_fail(SWIG_TypeError,"For director inputs, failed to convert input to " xName ".");
 }

%typemap(directorin, noblock=1, fragment="casadi_all") const xType& {
    if(!($input = swig::from< xType >($1))) %dirout_fail(SWIG_TypeError,"For director inputs, failed to convert input to " xName ".");
 }
 // Return-by-value; a target may replace these
#ifdef CASADI_TARGET_OUTPUT_TYPEMAPS
CASADI_TARGET_OUTPUT_TYPEMAPS(xName, xStubIn, xStubOut, xTsStub, xType)
#else
%typemap(out, doc=xName, pystub_out=xStubOut, tsstub_out=xTsStub, noblock=1, fragment="casadi_all") xType, const xType {
  if(!($result = swig::from< xType >($1))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output to type '" xName "'.");
}

// Return a const-ref behaves like return-by-value
%typemap(out, doc=xName, pystub_out=xStubOut, tsstub_out=xTsStub, noblock=1, fragment="casadi_all") const xType& {
  if(!($result = swig::from< xType >(*$1))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output to type '" xName "'.");
}
#endif

 // Enable dynamic dispatch
%typemap(typecheck, noblock=1, fragment="casadi_all") xType &OUTPUT {
  $1 = swig::check< xType >($input);
 }

// Alternative names
%apply xType &OUTPUT {xType &OUTPUT1};
%apply xType &OUTPUT {xType &OUTPUT2};
%apply xType &OUTPUT {xType &OUTPUT3};
%apply xType &OUTPUT {xType &OUTPUT4};
%apply xType &OUTPUT {xType &OUTPUT5};
%apply xType &OUTPUT {xType &OUTPUT6};

// Inputs marked INOUT are also returned by the function, ...
%typemap(argout,noblock=1,fragment="casadi_all") xType &INOUT {
  %append_output(swig::from(*$1));
 }

// ... but kept as inputs
%typemap(in, doc=xName, pystub_in=xStubIn, pystub_out=xStubOut, tsstub_in=xTsStub, tsstub_out=xTsStub, noblock=1, fragment="casadi_all") xType &INOUT (xType m) {
  $1 = &m;
  if (!swig::asref< xType >($input, &$1)) SWIG_exception_fail(SWIG_TypeError,"Failed to convert input to type '" xName "'.");
 }

 // ... also for dynamic dispatch
%typemap(typecheck, noblock=1, fragment="casadi_all") xType& INOUT {
  $1 = swig::check< xType >($input);
 }

// No arguments need to be freed
%typemap(freearg, noblock=1) xType& INOUT {}

// Alternative names
%apply xType &INOUT {xType &INOUT1};
%apply xType &INOUT {xType &INOUT2};
%apply xType &INOUT {xType &INOUT3};
%apply xType &INOUT {xType &INOUT4};
%apply xType &INOUT {xType &INOUT5};
%apply xType &INOUT {xType &INOUT6};
%enddef

 // All-in-one template instantiation + typemaps
%define %casadi_template(xName, xStubIn, xStubOut, xTsStub, xPrec, xType...)
%template() xType;
%casadi_input_typemaps(xName, xStubIn, xTsStub, xPrec, xType)
%casadi_output_typemaps(xName, xStubIn, xStubOut, xTsStub, %arg(xType))
%enddef

// std::vector<E> variant; xJSName names it where the target has vector classes
%define %casadi_template_vec(xName, xStubIn, xStubOut, xTsStub, xPrec, xJSName, xElemType...)
%casadi_template(xName, xStubIn, xStubOut, xTsStub, xPrec, std::vector< xElemType >)
#ifdef CASADI_TARGET_TEMPLATE_VEC
CASADI_TARGET_TEMPLATE_VEC(xJSName, xElemType)
#endif
%enddef

 // Input + output typemaps for an existing type.
%define %casadi_typemaps(xName, xStubIn, xStubOut, xTsStub, xPrec, xType...)
%casadi_input_typemaps(xName, xStubIn, xTsStub, xPrec, xType)
%casadi_output_typemaps(xName, xStubIn, xStubOut, xTsStub, xType)
%enddef


// Opti values as native numbers or arrays; targets without them override
%typemap(out, doc="double", pystub_out="float", noblock=1, fragment="casadi_all") casadi::native_DM {
  if(!($result = full_or_sparse($1, true))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output to type 'double'.");
}

// Order in typemap matching: Lower value means will be checked first

%define PREC_DICT 21 %enddef
%define PREC_SPARSITY 90 %enddef
%define PREC_IVector 92 %enddef
%define PREC_IVectorVector 92 %enddef
%define PREC_VECTOR 92 %enddef
%define PREC_SLICE 94 %enddef
%define PREC_IM 97 %enddef
%define PREC_DM 100 %enddef
%define PREC_DMVector 101 %enddef
%define PREC_DMVectorVector 101 %enddef
%define PREC_SX 103 %enddef
%define PREC_SXVector 103 %enddef
%define PREC_SXVectorVector 103 %enddef
%define PREC_MX 104 %enddef
%define PREC_MXVector 105 %enddef
%define PREC_MXVectorVector 106 %enddef
%define PREC_STRING 180 %enddef
%define PREC_FUNCTION 200 %enddef
%define PREC_GENERICTYPE 201 %enddef

#ifndef SWIGXML

// std::ostream & is redirected to casadi::uout()
%typemap(in, noblock=1, numinputs=0) std::ostream &stream {
  $1 = &casadi::uout();
}

// Add trailing newline in MATLAB and Octave

// Type names in error messages and docs; targets may define their own
#ifndef L_INT
#define L_INT "int"
#define L_BOOL "bool"
#define LPAIR(A,B) "(" A "," B ")"
  #define LL "["
  #define LR "]"
  #define L_DICT "dict"
  #define L_DOUBLE "float"
  #define LDICT(ARG) L_DICT ":" ARG
  #define L_STR "str"
#endif


#ifndef CASADI_TARGET_SCALAR_TYPEMAPS
%casadi_typemaps(L_STR, "builtins.str", "builtins.str", "string", PREC_STRING, std::string)
#endif
%casadi_template_vec(LL L_STR LR, "Sequence[builtins.str]", "list[builtins.str]", "string[]", PREC_VECTOR, StringVector, std::string)
%casadi_template(LL LL L_STR LR LR, "Sequence[Sequence[builtins.str]]", "list[list[builtins.str]]", "string[][]", PREC_VECTOR, std::vector<std::vector<std::string> >)

// PEP-484 stub aliases and preamble (-stubs)

#ifdef SWIG_STUBS_ENABLED

%insert("stubs_preamble") %{
# Overloads are ordered narrowest-first on purpose (`sin(float) -> float`
# before `sin(_DM) -> DM`), and DM/SX/MX compare elementwise; pyright
# reports both as overlap/override defects when this file itself is opened.
# test/python/pyright_stubs.py checks the stub with this line stripped.
# pyright: reportOverlappingOverload=false, reportIncompatibleMethodOverride=false
# `str` is also a method name on DM/SX/Opti/...; inside those class bodies a
# bare `str` annotation resolves to the method, so type tokens use builtins.str.
import builtins
import numpy as np
from numpy.typing import NDArray
from collections.abc import Callable
from typing import Protocol, runtime_checkable
_T = TypeVar("_T", DM, SX, MX)
# NZproxy is parameterised with Self, so its TypeVar is bounded, not
# value-constrained: a constrained one only accepts its exact constraints.
_TNZ = TypeVar("_TNZ", bound=DM | SX | MX)

@runtime_checkable
class _SupportsDM(Protocol):
    """Any object exposing a ``__DM__()`` method -- the conversion hook
    casadi uses to turn user-defined numeric types into a DM."""
    def __DM__(self) -> DM: ...

@runtime_checkable
class _SupportsTocsc(Protocol):
    """A scipy.sparse matrix or array, which casadi converts to a DM."""
    def tocsc(self) -> Any: ...

@runtime_checkable
class _SupportsSX(Protocol):
    def __SX__(self) -> SX: ...

@runtime_checkable
class _SupportsMX(Protocol):
    def __MX__(self) -> MX: ...

%}

// issue #2959: stub for casadi.ArrayInterface (a %pythoncode class)

// Matrix input aliases; SX and MX do not accept each other
%stub_alias_in(DM, bool | int | float | DM | Sparsity | _SupportsDM | _SupportsTocsc | Sequence[bool | int | float] | Sequence[Sequence[bool | int | float]] | NDArray[Any] | %arg(np.floating[Any]) | %arg(np.integer[Any]) | np.bool_)
%stub_alias_in(SX, SX | _SupportsSX | _DM)
%stub_alias_in(MX, MX | _SupportsMX | _DM)
// IM has no Python class: list what swig::asref<IM> accepts
%stub_alias_in(IM, int | Sequence[int] | Sequence[Sequence[int]] | NDArray[Any] | Sparsity | DM)
/* Out: swig::from(IM) hands back a DM (e.g. MX::mapping()). */
%stub_alias_out_key(IM, DM)
%stub_alias_in(SXElem, SXElem | float)

/* Bracketed doc tokens -- not valid SWIG identifiers. */
%stub_alias_in_key([float], Sequence[bool | int | float] | NDArray[Any])
%stub_alias_in_key([int],   Sequence[bool | int] | NDArray[Any])
%stub_alias_in_key([bool],  Sequence[bool | int] | NDArray[Any])

%stub_alias_in(Slice, Slice | int | slice)

// GenericType: precise input union, Any output
// Forward references stay unquoted inside | unions
%stub_alias_in(GenericType, bool | int | float | str | Function | Sequence[_GenericType] | %arg(Mapping[str, _GenericType]) | None)
%stub_alias_out_key(GenericType, Any)

// Matrix index; MX adds MX indices in its own %extend
%stub_alias_in(MIndex, int | slice | Slice | Sequence[bool | int] | NDArray[Any] | Sparsity | DM)

// CasadiMatrix operator overloads, narrowed per class (no MX + SX)
%define %stub_CasadiMatrix_binop(NAME)
%stub_overload_method_selftyped(NAME, DM, DM, other: _DM)
%stub_overload_method_selftyped(NAME, SX, SX, other: _SX)
%stub_overload_method_selftyped(NAME, MX, MX, other: _MX)
%enddef

%define %stub_CasadiMatrix_unop(NAME)
%stub_method0(NAME, Self)
%enddef

#else

%define %stub_CasadiMatrix_binop(NAME) %enddef
%define %stub_CasadiMatrix_unop(NAME)  %enddef

#endif

%casadi_typemaps("Sparsity", "Sparsity", "Sparsity", "Sparsity", PREC_SPARSITY, casadi::Sparsity)
%casadi_template_vec(LL "Sparsity" LR, "Sequence[Sparsity]", "list[Sparsity]", "Sparsity[]", PREC_SPARSITY, SparsityVector, casadi::Sparsity)
%casadi_template(LL LL "Sparsity"  LR  LR, "Sequence[Sequence[Sparsity]]", "list[list[Sparsity]]", "Sparsity[][]", PREC_SPARSITY, std::vector<std::vector< casadi::Sparsity> >)
%casadi_template(LDICT("Sparsity"), "Mapping[builtins.str, Sparsity]", "dict[builtins.str, Sparsity]", "Record<string, Sparsity>", PREC_SPARSITY, std::map<std::string, casadi::Sparsity >)
%casadi_template(LDICT(LL "Sparsity" LR), "Mapping[builtins.str, Sequence[Sparsity]]", "dict[builtins.str, list[Sparsity]]", "Record<string, Sparsity[]>", PREC_SPARSITY, std::map<std::string, std::vector<casadi::Sparsity > >)
%casadi_template(LPAIR(LDICT("Sparsity"),"[" L_STR "]"), "tuple[Mapping[builtins.str, Sparsity], Sequence[builtins.str]]", "tuple[dict[builtins.str, Sparsity], list[builtins.str]]", "[Record<string, Sparsity>, string[]]", PREC_SPARSITY, std::pair<std::map<std::string, casadi::Sparsity >, std::vector<std::string> >)
#ifndef CASADI_TARGET_SCALAR_TYPEMAPS
%casadi_typemaps(L_BOOL, "bool", "bool", "boolean", SWIG_TYPECHECK_BOOL, bool)
%casadi_template("[" L_BOOL "]", "Sequence[bool | int] | NDArray[Any]", "list[bool]", "boolean[]", SWIG_TYPECHECK_BOOL, std::vector<bool>)
%casadi_template("[[" L_BOOL "]]", "Sequence[Sequence[bool | int]] | Sequence[NDArray[Any]]", "list[list[bool]]", "boolean[][]", SWIG_TYPECHECK_BOOL, std::vector<std::vector<bool> >)
%casadi_typemaps( L_INT , "int", "int", "bigint", SWIG_TYPECHECK_INTEGER, casadi_int)
#endif

#ifndef LABEL
#define LABEL LPAIR("int","int")
#endif
%casadi_template(LABEL, "tuple[int, int]", "tuple[int, int]", "[bigint, bigint]", SWIG_TYPECHECK_INTEGER, std::pair<casadi_int,casadi_int>)
#undef LABEL
%casadi_template_vec("[" L_INT "]", "Sequence[bool | int] | NDArray[Any]", "list[int]", "bigint[]", PREC_IVector, IntVector, casadi_int)
%casadi_template(LL "[" L_INT "]" LR, "Sequence[Sequence[bool | int] | NDArray[Any]]", "list[list[int]]", "bigint[][]", PREC_IVectorVector, std::vector<std::vector<casadi_int> >)
#ifndef CASADI_TARGET_SCALAR_TYPEMAPS
%casadi_typemaps(L_DOUBLE, "float", "float", "number", SWIG_TYPECHECK_DOUBLE, double)
#endif
%casadi_template_vec("[" L_DOUBLE "]", "Sequence[bool | int | float] | NDArray[Any]", "list[float]", "number[]", SWIG_TYPECHECK_DOUBLE, DoubleVector, double)
%casadi_template(LL "[" L_DOUBLE "]" LR, "Sequence[Sequence[bool | int | float] | NDArray[Any]]", "list[list[float]]", "number[][]", SWIG_TYPECHECK_DOUBLE, std::vector<std::vector<double> >)
// Most specific first: DM, SX, MX (d.ts overload order)
%casadi_typemaps("DM", "DM", "DM", "DM", PREC_DM, casadi::Matrix<double>)
%casadi_template_vec(LL "DM" LR, "Sequence[DM]", "list[DM]", "DM[]", PREC_DMVector, DMVector, casadi::Matrix<double>)
%casadi_template(LL LL "DM" LR LR, "Sequence[Sequence[DM]]", "list[list[DM]]", "DM[][]", PREC_DMVectorVector, std::vector<std::vector< casadi::Matrix<double> > >)
%casadi_template(LDICT("DM"), "Mapping[builtins.str, DM]", "dict[builtins.str, DM]", "Record<string, DM>", PREC_DM, std::map<std::string, casadi::Matrix<double> >)
%casadi_typemaps("SXElem", "SXElem", "SXElem", "SXElem", PREC_SX, casadi::SXElem)
%casadi_template(LL "SXElem" LR, "Sequence[SXElem]", "list[SXElem]", "SXElem[]", PREC_SXVector, std::vector<casadi::SXElem>)
%casadi_typemaps("SX", "SX", "SX", "SX", PREC_SX, casadi::Matrix<casadi::SXElem>)
%casadi_template_vec(LL "SX" LR, "Sequence[SX]", "list[SX]", "SX[]", PREC_SXVector, SXVector, casadi::Matrix<casadi::SXElem>)
%casadi_template(LL LL "SX" LR LR, "Sequence[Sequence[SX]]", "list[list[SX]]", "SX[][]", PREC_SXVectorVector, std::vector<std::vector< casadi::Matrix<casadi::SXElem> > >)
%casadi_template(LDICT("SX"), "Mapping[builtins.str, SX]", "dict[builtins.str, SX]", "Record<string, SX>", PREC_SX, std::map<std::string, casadi::Matrix<casadi::SXElem> >)
%casadi_typemaps("MX", "MX", "MX", "MX", PREC_MX, casadi::MX)
%casadi_template_vec(LL "MX" LR, "Sequence[MX]", "list[MX]", "MX[]", PREC_MXVector, MXVector, casadi::MX)
%casadi_template(LL LL "MX" LR LR, "Sequence[Sequence[MX]]", "list[list[MX]]", "MX[][]", PREC_MXVectorVector, std::vector<std::vector<casadi::MX> >)
%casadi_template(LDICT("MX"), "Mapping[builtins.str, MX]", "dict[builtins.str, MX]", "Record<string, MX>", PREC_MX, std::map<std::string, casadi::MX>)
%casadi_template(LPAIR("MX","MX"), "tuple[MX, MX]", "tuple[MX, MX]", "[MX, MX]", PREC_MXVector, std::pair<casadi::MX, casadi::MX>)
%casadi_typemaps("IM", "IM", "IM", "IM", PREC_IM, casadi::Matrix<casadi_int>)
%casadi_typemaps("GenericType", "GenericType", "GenericType", "GenericType", PREC_GENERICTYPE, casadi::GenericType)
%casadi_template(LL "GenericType" LR, "Sequence[GenericType]", "list[GenericType]", "GenericType[]", PREC_GENERICTYPE, std::vector<casadi::GenericType>)
%casadi_template(LL LL "GenericType" LR LR, "Sequence[Sequence[GenericType]]", "list[list[GenericType]]", "GenericType[][]", PREC_GENERICTYPE, std::vector<std::vector<casadi::GenericType> >)
%casadi_typemaps("Slice", "Slice", "Slice", "Slice", PREC_SLICE, casadi::Slice)
%casadi_typemaps("Function", "Function", "Function", "Function", PREC_FUNCTION, casadi::Function)
%casadi_template_vec(LL "Function" LR, "Sequence[Function]", "list[Function]", "Function[]", PREC_FUNCTION, FunctionVector, casadi::Function)
%casadi_template(LPAIR("Function","Function"), "tuple[Function, Function]", "tuple[Function, Function]", "[Function, Function]", PREC_FUNCTION, std::pair<casadi::Function, casadi::Function>)
%casadi_template(L_DICT, "Mapping[builtins.str, GenericType]", "dict[builtins.str, GenericType]", "Record<string, GenericType>", PREC_DICT, std::map<std::string, casadi::GenericType>)
%casadi_template(LDICT(LL L_STR LR), "Mapping[builtins.str, Sequence[builtins.str]]", "dict[builtins.str, list[builtins.str]]", "Record<string, string[]>", PREC_DICT, std::map<std::string, std::vector<std::string> >)


#undef L_INT
#undef L_BOOL
#undef LPAIR
#undef L_DOUBLE
#undef L_DICT
#undef LL
#undef LR
#undef L_STR

// TypeScript aliases for casadi typedefs used in signatures

// Matlab is index-1 based

#endif // SWIGXML


// Init hooks


// Strip leading casadi_ unless followed by ML/int
%rename("%(regex:/casadi_(?!ML|int\\b)(.*)/\\1/)s") "";

%rename(row) get_row;
%rename(colind) get_colind;
%rename(sparsity) get_sparsity;
%rename(nonzeros) get_nonzeros;
%rename(elements) get_elements;

// Explicit conversion to double and casadi_int
%rename(to_double) operator double;
%rename(to_int) operator casadi_int;
%rename(to_DM) operator Matrix<double>;


// Indexing helpers of the matrix classes; targets define their own
#ifndef %matrix_helpers
%define %matrix_helpers(Type, xIndex) %enddef
#endif


#ifdef SWIGPYTHON
%include "python/casadi_python_decl.i"
#elif defined(SWIGMATLAB)
%include "matlab/casadi_matlab_decl.i"
#elif defined(SWIGJULIA)
%include "julia/casadi_julia_decl.i"
#elif defined(SWIGWASMJS)
%include "wasm-js/casadi_wasmjs_decl.i"
#elif defined(SWIGJAVA)
%include "java/casadi_java_decl.i"
#elif defined(SWIGRUST)
%include "rust/casadi_rust_decl.i"
#elif defined(SWIGDOTNET)
%include "dotnet/casadi_dotnet_decl.i"
#endif
%include <casadi/core/printable.hpp>

%include <casadi/core/generic_shared.hpp>

%template(GenSharedObject) casadi::GenericShared<casadi::SharedObject, casadi::SharedObjectInternal>;
%template(GenWeakRef) casadi::GenericWeakRef<casadi::SharedObject, casadi::SharedObjectInternal>;

%include <casadi/core/shared_object.hpp>
%include <casadi/core/casadi_misc.hpp>
%include <casadi/core/casadi_common.hpp>
%include <casadi/core/generic_type.hpp>
%include <casadi/core/calculus.hpp>
%include <casadi/core/sparsity_interface.hpp>
%include <casadi/core/sparsity.hpp>

// Logic for pickling

%include <casadi/core/slice.hpp>


%include <casadi/core/generic_matrix.hpp>

%template(GenDM)        casadi::GenericMatrix<casadi::Matrix<double> >;
%template(GenSX)             casadi::GenericMatrix<casadi::Matrix<casadi::SXElem> >;
%template(GenMX)             casadi::GenericMatrix<casadi::MX>;

%include <casadi/core/generic_expression.hpp>

// Interface-internal flags, not module constants
%ignore IS_GLOBAL; %ignore IS_MEMBER; %ignore IS_SPARSITY; %ignore IS_DMATRIX;
%ignore IS_IMATRIX; %ignore IS_SX; %ignore IS_MX; %ignore IS_DOUBLE; %ignore FLAG;
#define IS_GLOBAL   0x1
#define IS_MEMBER   0x10
#define IS_SPARSITY 0x100
#define IS_DMATRIX  0x1000
#define IS_IMATRIX  0x10000
#define IS_SX       0x100000
#define IS_MX       0x1000000
#define IS_DOUBLE   0x10000000

%define SPARSITY_INTERFACE_FUN_BASE(DECL, FLAG, M)
#if FLAG & IS_MEMBER

 DECL M casadi_horzcat(const std::vector< M > &v) {
  return horzcat(v);
 }
 DECL M casadi_vertcat(const std::vector< M > &v) {
 return vertcat(v);
 }
 DECL std::vector< M >
 casadi_horzsplit(const M& v, const std::vector<casadi_int>& offset) {
 return horzsplit(v, offset);
 }
 DECL std::vector< M > casadi_horzsplit(const M& v, casadi_int incr=1) {
 return horzsplit(v, incr);
 }
 DECL std::vector< M > casadi_horzsplit_n(const M& v, casadi_int n) {
 return horzsplit_n(v, n);
 }
 DECL std::vector< M >
 casadi_vertsplit(const M& v, const std::vector<casadi_int>& offset) {
 return vertsplit(v, offset);
 }
 DECL std::vector<casadi_int >
 casadi_offset(const std::vector< M > &v, bool vert=true) {
 return offset(v, vert);
 }
 DECL std::vector< M >
 casadi_vertsplit(const M& v, casadi_int incr=1) {
 return vertsplit(v, incr);
 }
 DECL std::vector< M >
 casadi_vertsplit_n(const M& v, casadi_int n) {
 return vertsplit_n(v, n);
 }
 DECL M casadi_blockcat(const M& A, const M& B, const M& C, const M& D) {
 return vertcat(horzcat(A, B), horzcat(C, D));
 }
 DECL std::vector< std::vector< M > >
 casadi_blocksplit(const M& x, const std::vector<casadi_int>& vert_offset,
 const std::vector<casadi_int>& horz_offset) {
 return blocksplit(x, vert_offset, horz_offset);
 }
 DECL std::vector< std::vector< M > >
 casadi_blocksplit(const M& x, casadi_int vert_incr=1, casadi_int horz_incr=1) {
 return blocksplit(x, vert_incr, horz_incr);
 }
 DECL M casadi_diagcat(const std::vector< M > &A) {
 return diagcat(A);
 }
 DECL std::vector< M >
 casadi_diagsplit(const M& x, const std::vector<casadi_int>& output_offset1,
 const std::vector<casadi_int>& output_offset2) {
 return diagsplit(x, output_offset1, output_offset2);
 }
 DECL std::vector< M >
 casadi_diagsplit(const M& x, const std::vector<casadi_int>& output_offset) {
 return diagsplit(x, output_offset);
 }
 DECL std::vector< M > casadi_diagsplit(const M& x, casadi_int incr=1) {
 return diagsplit(x, incr);
 }
 DECL std::vector< M >
 casadi_diagsplit(const M& x, casadi_int incr1, casadi_int incr2) {
 return diagsplit(x, incr1, incr2);
 }
 DECL M casadi_veccat(const std::vector< M >& x) {
 return veccat(x);
 }
 DECL M casadi_mtimes(const M& x, const M& y,
                      const std::string& blas = "reference") {
 return mtimes(x, y, blas);
 }
 DECL M casadi_mtimes(const std::vector< M > &args) {
 return mtimes(args);
 }
 DECL M casadi_mtimes(const std::vector< M > &args, const std::string& blas) {
 return mtimes(args, blas);
 }
 DECL M casadi_mac(const M& X, const M& Y, const M& Z,
                   const std::string& blas = "reference") {
 return mac(X, Y, Z, blas);
 }
 DECL M casadi_transpose(const M& X) {
 return X.T();
 }
 DECL M casadi_vec(const M& a) {
 return vec(a);
 }
 DECL M casadi_reshape(const M& a, casadi_int nrow, casadi_int ncol) {
 return reshape(a, nrow, ncol);
 }
 DECL M casadi_reshape(const M& a, std::pair<casadi_int, casadi_int> rc) {
 return reshape(a, rc.first, rc.second);
 }
 DECL M casadi_reshape(const M& a, const Sparsity& sp) {
 return reshape(a, sp);
 }
 DECL M casadi_sparsity_cast(const M& a, const Sparsity& sp) {
 return sparsity_cast(a, sp);
 }
 DECL casadi_int casadi_sprank(const M& A) {
 return sprank(A);
 }
 DECL casadi_int casadi_norm_0_mul(const M& x, const M& y) {
 return norm_0_mul(x, y);
 }
 DECL M casadi_triu(const M& a, bool includeDiagonal=true) {
 return triu(a, includeDiagonal);
 }
 DECL M casadi_tril(const M& a, bool includeDiagonal=true) {
 return tril(a, includeDiagonal);
 }
 DECL M casadi_kron(const M& a, const M& b) {
 return kron(a, b);
 }
 DECL M casadi_repmat(const M& A, casadi_int n, casadi_int m=1) {
 return repmat(A, n, m);
 }
 DECL M casadi_repmat(const M& A, const std::pair<casadi_int, casadi_int>& rc) {
 return repmat(A, rc.first, rc.second);
 }
 DECL M casadi_sum2(const M& x) {
 return sum2(x);
 }
 DECL M casadi_sum1(const M& x) {
 return sum1(x);
 }
#endif
#ifdef CASADI_TARGET_SPARSITY_INTERFACE_FUN
CASADI_TARGET_SPARSITY_INTERFACE_FUN(DECL, FLAG, M)
#endif
%enddef

%define SPARSITY_INTERFACE_ALL(DECL, FLAG)
SPARSITY_INTERFACE_FUN(DECL, (FLAG | IS_SPARSITY), Sparsity)
SPARSITY_INTERFACE_FUN(DECL, (FLAG | IS_MX), MX)
SPARSITY_INTERFACE_FUN(DECL, (FLAG | IS_DMATRIX), Matrix<double>)
SPARSITY_INTERFACE_FUN(DECL, (FLAG | IS_SX), Matrix<SXElem>)
%enddef

  %define SPARSITY_INTERFACE_FUN(DECL, FLAG, M)
    SPARSITY_INTERFACE_FUN_BASE(DECL, FLAG, M)
    #if FLAG & IS_MEMBER
#ifdef CASADI_TARGET_SUM_FUN
CASADI_TARGET_SUM_FUN(DECL, M)
#else
      DECL M casadi_sum(const M& x, casadi_int dim) {
        if (dim==0) return sum1(x);
        if (dim==1) return sum2(x);
        casadi_error(
          "Expected sum(A,0), sum(A,1), sum(A,\"all\") got " + casadi::str(dim) + " instead.");
      }
      DECL M casadi_sum(const M& x) {
        return sum(x);
      }
#endif
    #endif
  %enddef

%define GENERIC_MATRIX_FUN(DECL, FLAG, M)
#if FLAG & IS_MEMBER
DECL M casadi_mpower(const M& x, const M& n) {
  return mpower(x, n);
}

DECL M casadi_mrdivide(const M& x, const M& y) {
  return mrdivide(x, y);
}

DECL M casadi_mldivide(const M& x, const M& y) {
  return mldivide(x, y);
}

DECL std::vector< M > casadi_symvar(const M& x) {
  return symvar(x);
}

DECL M casadi_transform(const M& x, const casadi::Dict& opts = casadi::Dict()) {
  return transform(x, opts);
}

DECL M casadi_transform(const M& x, const std::vector< std::vector< casadi::GenericType > >& passes, const casadi::Dict& opts = casadi::Dict()) {
  return transform(x, passes, opts);
}

DECL std::vector< M > casadi_transform(const std::vector< M >& x, const casadi::Dict& opts = casadi::Dict()) {
  return transform(x, opts);
}

DECL std::vector< M > casadi_transform(const std::vector< M >& x, const std::vector< std::vector< casadi::GenericType > >& passes, const casadi::Dict& opts = casadi::Dict()) {
  return transform(x, passes, opts);
}

DECL M casadi_bilin(const M& A, const M& x, const M& y) {
  return bilin(A, x, y);
}

DECL M casadi_bilin(const M& A, const M& x) {
  return bilin(A, x);
}

DECL M casadi_rank1(const M& A, const M& alpha, const M& x, const M& y) {
  return rank1(A, alpha, x, y);
}

DECL M casadi_sumsqr(const M& X) {
  return sumsqr(X);
}

DECL M casadi_linspace(const M& a, const M& b, casadi_int nsteps) {
  return linspace(a, b, nsteps);
}

DECL M casadi_logsumexp(const M& a) {
  return logsumexp(a);
}

DECL M casadi_logsumexp(const M& a, const M& margin) {
  return logsumexp(a, margin);
}

DECL M casadi_interp1d(const std::vector<double>& x, const M&v,
        const std::vector<double>& xq, const std::string& mode="linear", bool equidistant=false) {
  return interp1d(x, v, xq, mode, equidistant);
}

DECL M casadi_soc(const M& x, const M& y) {
  return soc(x, y);
}

DECL M casadi_cross(const M& a, const M& b, casadi_int dim = -1) {
  return cross(a, b, dim);
}

DECL M casadi_skew(const M& a) {
  return skew(a);
}

DECL M casadi_inv_skew(const M& a) {
  return inv_skew(a);
}

DECL M casadi_det(const M& A) {
  return det(A);
}

DECL M casadi_det(const M& A, const std::string& lsolver,
                      const casadi::Dict& opts = casadi::Dict()) {
  return det(A, lsolver, opts);
}

DECL M casadi_inv_minor(const M& A) {
  return inv_minor(A);
}

DECL M casadi_inv(const M& A) {
  return inv(A);
}

DECL M casadi_inv(const M& A, const std::string& lsolver,
                      const casadi::Dict& opts = casadi::Dict()) {
  return inv(A, lsolver, opts);
}

DECL M casadi_trace(const M& a) {
  return trace(a);
}

DECL M casadi_tril2symm(const M& a) {
  return tril2symm(a);
}

DECL M casadi_triu2symm(const M& a) {
  return triu2symm(a);
}

DECL M casadi_norm_fro(const M& x) {
  return norm_fro(x);
}

DECL M casadi_norm_2(const M& x) {
  return norm_2(x);
}

DECL M casadi_norm_1(const M& x) {
  return norm_1(x);
}

DECL M casadi_norm_inf(const M& x) {
  return norm_inf(x);
}

DECL M casadi_dot(const M& x, const M& y) {
  return dot(x, y);
}

DECL M casadi_nullspace(const M& A) {
  return nullspace(A);
}

DECL M casadi_polyval(const M& p, const M& x) {
  return polyval(p, x);
}

DECL M casadi_diag(const M& A) {
  return diag(A);
}

DECL M casadi_unite(const M& A, const M& B) {
  return unite(A, B);
}

DECL M casadi_densify(const M& x) {
  return densify(x);
}

DECL M casadi_densify(const M& x, const M& val) {
  return densify(x, val);
}

DECL M casadi_project(const M& A, const Sparsity& sp, bool intersect=false) {
  return project(A, sp, intersect);
}

DECL M casadi_if_else(const M& cond, const M& if_true,
                    const M& if_false, bool short_circuit=false) {
  return if_else(cond, if_true, if_false, short_circuit);
}

DECL M casadi_conditional(const M& ind, const std::vector< M > &x,
                        const M& x_default, bool short_circuit=false) {
  return conditional(ind, x, x_default, short_circuit);
}

DECL bool casadi_depends_on(const M& f, const M& arg) {
  return depends_on(f, arg);
}

DECL bool casadi_contains(const std::vector<M>& v, const M& n) {
  return contains(v, n);
}

DECL bool casadi_contains_all(const std::vector<M>& v, const std::vector<M>& n) {
  return contains_all(v, n);
}

DECL bool casadi_contains_any(const std::vector<M>& v, const std::vector<M>& n) {
  return contains_any(v, n);
}

DECL M casadi_solve(const M& A, const M& b) {
  return solve(A, b);
}

DECL M casadi_solve(const M& A, const M& b,
                       const std::string& lsolver,
                       const casadi::Dict& opts = casadi::Dict()) {
  return solve(A, b, lsolver, opts);
}

DECL M casadi_pinv(const M& A) {
  return pinv(A);
}

DECL M casadi_pinv(const M& A, const std::string& lsolver,
                      const casadi::Dict& opts = casadi::Dict()) {
  return pinv(A, lsolver, opts);
}

DECL M casadi_expm_const(const M& A, const M& t) {
  return expm_const(A, t);
}

DECL M casadi_expm(const M& A) {
  return expm(A);
}

DECL M casadi_jacobian(const M &ex, const M &arg, const Dict& opts=Dict()) {
  return jacobian(ex, arg, opts);
}

DECL M casadi_jtimes(const M& ex, const M& arg, const M& v, bool tr=false, const Dict& opts=Dict()) {
  return jtimes(ex, arg, v, tr);
}

DECL M casadi_linearize(const M& f, const M& x, const M& x0) {
  return linearize(f, x, x0);
}

DECL std::vector<bool> casadi_which_depends(const M& expr, const M& var,
                                            casadi_int order=1, bool tr=false) {
  return which_depends(expr, var, order, tr);
}

DECL Sparsity casadi_jacobian_sparsity(const M& f, const M& x) {
  return jacobian_sparsity(f, x);
}

DECL bool casadi_is_linear(const M& expr, const M& var) {
  return is_linear(expr, var);
}

DECL bool casadi_is_quadratic(const M& expr, const M& var) {
  return is_quadratic(expr, var);
}

DECL M casadi_gradient(const M &ex, const M &arg, const Dict& opts=Dict()) {
  return gradient(ex, arg, opts);
}

DECL M casadi_tangent(const M &ex, const M &arg, const Dict& opts=Dict()) {
  return tangent(ex, arg, opts);
}

DECL M casadi_hessian(const M& ex, const M& arg, M& OUTPUT1, const casadi::Dict& opts = casadi::Dict()) {
  return hessian(ex, arg, OUTPUT1, opts);
}

DECL void casadi_quadratic_coeff(const M& ex, const M& arg, M& OUTPUT1, M& OUTPUT2, M& OUTPUT3, bool check=true) {
  quadratic_coeff(ex, arg, OUTPUT1, OUTPUT2, OUTPUT3, check);
}

DECL void casadi_linear_coeff(const M& ex, const M& arg, M& OUTPUT1, M& OUTPUT2, bool check=true) {
  linear_coeff(ex, arg, OUTPUT1, OUTPUT2, check);
}

DECL casadi_int casadi_n_nodes(const M& A) {
  return n_nodes(A);
}

DECL std::string casadi_print_operator(const M& xb,
                                                  const std::vector<std::string>& args) {
  return print_operator(xb, args);
}
DECL M casadi_repsum(const M& A, casadi_int n, casadi_int m=1) {
  return repsum(A, n, m);
}
DECL M casadi_diff(const M& A, casadi_int n=1, casadi_index axis=-1) {
  return diff(A, n, axis);
}
DECL M casadi_cumsum(const M& A, casadi_index axis=-1) {
  return cumsum(A, axis);
}
DECL M casadi_einstein(const M& A, const M& B, const M& C,
  const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b, const std::vector<casadi_int>& dim_c,
  const std::vector<casadi_int>& a, const std::vector<casadi_int>& b, const std::vector<casadi_int>& c) {
  return einstein(A, B, C, dim_a, dim_b, dim_c, a, b, c);
}
DECL M casadi_einstein(const M& A, const M& B,
  const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b, const std::vector<casadi_int>& dim_c,
  const std::vector<casadi_int>& a, const std::vector<casadi_int>& b, const std::vector<casadi_int>& c) {
  return einstein(A, B, dim_a, dim_b, dim_c, a, b, c);
}
DECL M casadi_mmin(const M& x) { return mmin(x); }
DECL M casadi_mmax(const M& x) { return mmax(x); }
DECL casadi::DM casadi_evalf(const M& x) {
  return evalf(x);
}
DECL void casadi_separate_linear(const M &expr,
      const M &sym_lin, const M &sym_const,
      M& OUTPUT1, M& OUTPUT2, M& OUTPUT3) {
  separate_linear(expr, sym_lin, sym_const, OUTPUT1, OUTPUT2, OUTPUT3);
}
DECL void casadi_separate_linear(const M &expr,
  const std::vector<M> &sym_lin, const std::vector<M> &sym_const,
  M& OUTPUT1, M& OUTPUT2, M& OUTPUT3) {
separate_linear(expr, sym_lin, sym_const, OUTPUT1, OUTPUT2, OUTPUT3);
}
#endif // FLAG & IS_MEMBER

#if FLAG & IS_GLOBAL
DECL std::vector<M> casadi_cse(const std::vector<M>& e) {
  return cse(e);
}
DECL M casadi_cse(const M& e) {
  return cse(e);
}

DECL void casadi_extract_parametric(const M &expr, const M& par,
        M& OUTPUT1, std::vector<M>& OUTPUT2, std::vector<M>& OUTPUT3, const Dict& opts=Dict()) {
  extract_parametric(expr, par, OUTPUT1, OUTPUT2, OUTPUT3, opts);
}
DECL void casadi_extract_parametric(const M &expr, const std::vector<M>& par,
        M& OUTPUT1, std::vector<M>& OUTPUT2, std::vector<M>& OUTPUT3, const Dict& opts=Dict()) {
  extract_parametric(expr, par, OUTPUT1, OUTPUT2, OUTPUT3, opts);
}
DECL void casadi_extract_parametric(const std::vector<M> &expr, const M& par,
        std::vector<M>& OUTPUT1, std::vector<M>& OUTPUT2, std::vector<M>& OUTPUT3, const Dict& opts=Dict()) {
  extract_parametric(expr, par, OUTPUT1, OUTPUT2, OUTPUT3, opts);
}
DECL void casadi_extract_parametric(const std::vector<M> &expr, const std::vector<M>& par,
        std::vector<M>& OUTPUT1, std::vector<M>& OUTPUT2, std::vector<M>& OUTPUT3, const Dict& opts=Dict()) {
  extract_parametric(expr, par, OUTPUT1, OUTPUT2, OUTPUT3, opts);
}

DECL std::vector<std::vector< M > >
casadi_forward(const std::vector< M > &ex, const std::vector< M > &arg,
               const std::vector<std::vector< M > > &v,
               const Dict& opts = Dict()) {
  return forward(ex, arg, v, opts);
}

DECL std::vector<std::vector< M > >
casadi_reverse(const std::vector< M > &ex, const std::vector< M > &arg,
               const std::vector<std::vector< M > > &v,
               const Dict& opts = Dict()) {
  return reverse(ex, arg, v, opts);
}

#if FLAG & (IS_SX | IS_MX)
DECL std::string casadi_export_graph(const std::vector<M>& expressions, const Dict& opts=Dict()) {
  return export_graph(expressions, opts);
}

DECL std::string casadi_export_graph(const M& expression, const Dict& opts=Dict()) {
  return export_graph(expression, opts);
}

DECL void casadi_export_graph(const std::vector<M>& expressions, const std::string& fname,
                             const Dict& opts=Dict()) {
  export_graph(expressions, fname, opts);
}

DECL void casadi_export_graph(const M& expression, const std::string& fname,
                             const Dict& opts=Dict()) {
  export_graph(expression, fname, opts);
}
#endif

DECL M casadi_substitute(const M& ex, const M& v, const M& vdef) {
  return substitute(ex, v, vdef);
}

DECL std::vector< M > casadi_substitute(const std::vector< M >& ex,
                                         const std::vector< M >& v,
                                         const std::vector< M >& vdef) {
  return substitute(ex, v, vdef);
}

DECL void casadi_substitute_inplace(const std::vector< M >& v,
                                      std::vector< M >& INOUT1,
                                      std::vector< M >& INOUT2,
                                      bool reverse=false) {
  return substitute_inplace(v, INOUT1, INOUT2, reverse);
}

DECL void casadi_extract(const std::vector< M >& ex,
    std::vector< M >& OUTPUT1,
    std::vector< M >& OUTPUT2,
    std::vector< M >& OUTPUT3,
    const Dict& opts = Dict()) {
  OUTPUT1 = ex;
  extract(OUTPUT1, OUTPUT2, OUTPUT3, opts);
}

DECL void casadi_shared(const std::vector< M >& ex,
                               std::vector< M >& OUTPUT1,
                               std::vector< M >& OUTPUT2,
                               std::vector< M >& OUTPUT3,
                               const std::string& v_prefix="v_",
                               const std::string& v_suffix="") {
  OUTPUT1 = ex;
  shared(OUTPUT1, OUTPUT2, OUTPUT3, v_prefix, v_suffix);
}

DECL M casadi_blockcat(const std::vector< std::vector< M > > &v) {
 return blockcat(v);
}
#endif // FLAG & IS_GLOBAL
%enddef

%define GENERIC_MATRIX_ALL(DECL, FLAG)
GENERIC_MATRIX_FUN(DECL, (FLAG | IS_MX), MX)
GENERIC_MATRIX_FUN(DECL, (FLAG | IS_DMATRIX), Matrix<double>)
GENERIC_MATRIX_FUN(DECL, (FLAG | IS_SX), Matrix<SXElem>)
%enddef

%define GENERIC_EXPRESSION_FUN(DECL, FLAG, M)
#if FLAG & IS_MEMBER
DECL M casadi_plus(const M& x, const M& y) { return x+y; }
DECL M casadi_minus(const M& x, const M& y) { return x-y; }
DECL M casadi_times(const M& x, const M& y) { return x*y; }
DECL M casadi_rdivide(const M& x, const M& y) { return x/y; }
DECL M casadi_ldivide(const M& x, const M& y) { return y/x; }
DECL M casadi_lt(const M& x, const M& y) { return x<y; }
DECL M casadi_le(const M& x, const M& y) { return x<=y; }
DECL M casadi_gt(const M& x, const M& y) { return x>y; }
DECL M casadi_ge(const M& x, const M& y) { return x>=y; }
DECL M casadi_eq(const M& x, const M& y) { return x==y; }
DECL M casadi_ne(const M& x, const M& y) { return x!=y; }
DECL M casadi_and(const M& x, const M& y) { return x&&y; }
DECL M casadi_or(const M& x, const M& y) { return x||y; }
DECL M casadi_not(const M& x) { return !x; }
DECL M casadi_abs(const M& x) { return fabs(x); }
DECL M casadi_sqrt(const M& x) { return sqrt(x); }
DECL M casadi_sin(const M& x) { return sin(x); }
DECL M casadi_cos(const M& x) { return cos(x); }
DECL M casadi_tan(const M& x) { return tan(x); }
DECL M casadi_atan(const M& x) { return atan(x); }
DECL M casadi_asin(const M& x) { return asin(x); }
DECL M casadi_acos(const M& x) { return acos(x); }
DECL M casadi_tanh(const M& x) { return tanh(x); }
DECL M casadi_sinh(const M& x) { return sinh(x); }
DECL M casadi_cosh(const M& x) { return cosh(x); }
DECL M casadi_atanh(const M& x) { return atanh(x); }
DECL M casadi_asinh(const M& x) { return asinh(x); }
DECL M casadi_acosh(const M& x) { return acosh(x); }
DECL M casadi_exp(const M& x) { return exp(x); }
DECL M casadi_log(const M& x) { return log(x); }
DECL M casadi_log10(const M& x) { return log10(x); }
DECL M casadi_log1p(const M& x) { return log1p(x); }
DECL M casadi_expm1(const M& x) { return expm1(x); }
DECL M casadi_floor(const M& x) { return floor(x); }
DECL M casadi_ceil(const M& x) { return ceil(x); }
DECL M casadi_erf(const M& x) { return erf(x); }
DECL M casadi_erfinv(const M& x) { using casadi::erfinv; return erfinv(x); }
DECL M casadi_sign(const M& x) { using casadi::sign; return sign(x); }
DECL M casadi_power(const M& x, const M& n) { return pow(x, n); }
DECL M casadi_mod(const M& x, const M& y) { return fmod(x, y); }
DECL M casadi_fmod(const M& x, const M& y) { return fmod(x, y); }
DECL M casadi_remainder(const M& x, const M& y) { return remainder(x, y); }
DECL M casadi_atan2(const M& x, const M& y) { return atan2(x, y); }
DECL M casadi_fmin(const M& x, const M& y) { return fmin(x, y); }
DECL M casadi_fmax(const M& x, const M& y) { return fmax(x, y); }
DECL M casadi_hypot(const M& x, const M& y) { return hypot(x, y); }
DECL M casadi_simplify(const M& x) { using casadi::simplify; return simplify(x); }
DECL bool casadi_is_equal(const M& x, const M& y, casadi_int depth=0) { using casadi::is_equal; return is_equal(x, y, depth); }
DECL M casadi_copysign(const M& x, const M& y) { return copysign(x, y); }
DECL M casadi_constpow(const M& x, const M& y) { using casadi::constpow; return constpow(x, y); }
#endif // FLAG & IS_MEMBER
%enddef

%define GENERIC_EXPRESSION_ALL(DECL, FLAG)
GENERIC_EXPRESSION_FUN(DECL, (FLAG | IS_MX), MX)
GENERIC_EXPRESSION_FUN(DECL, (FLAG | IS_DMATRIX), Matrix<double>)
GENERIC_EXPRESSION_FUN(DECL, (FLAG | IS_SX), Matrix<SXElem>)
GENERIC_EXPRESSION_FUN(DECL, (FLAG | IS_DOUBLE), double)
%enddef

%define MATRIX_FUN(DECL, FLAG, M)
#if FLAG & IS_MEMBER
DECL M casadi_all(const M& x) {
  return all(x);
}

DECL M casadi_any(const M& x) {
  return any(x);
}

DECL M casadi_adj(const M& A) {
  return adj(A);
}

DECL M casadi_minor(const M& x, casadi_int i, casadi_int j) {
  return minor(x, i, j);
}

DECL M casadi_cofactor(const M& x, casadi_int i, casadi_int j) {
  return cofactor(x, i, j);
}

DECL void casadi_qr(const M& A, M& OUTPUT1, M& OUTPUT2) {
  return qr(A, OUTPUT1, OUTPUT2);
}

DECL void casadi_qr_sparse(const M& A, M& OUTPUT1, M& OUTPUT2, M& OUTPUT3,
          std::vector<casadi_int>& OUTPUT4, std::vector<casadi_int>& OUTPUT5, bool amd=true) {
  return qr_sparse(A, OUTPUT1, OUTPUT2, OUTPUT3, OUTPUT4, OUTPUT5, amd);
}

DECL M casadi_qr_solve(const M& b, const M& v, const M& r, const M& beta,
                       const std::vector<casadi_int>& prinv,
                       const std::vector<casadi_int>& pc, bool tr=false) {
  return qr_solve(b, v, r, beta, prinv, pc, tr);
}

DECL void casadi_ldl(const M& A, M& OUTPUT1, M& OUTPUT2, std::vector<casadi_int>& OUTPUT3, bool amd=true) {
  return ldl(A, OUTPUT1, OUTPUT2, OUTPUT3, amd);
}

DECL M casadi_ldl_solve(const M& b, const M& D, const M& LT, const std::vector<casadi_int>& p) {
  return ldl_solve(b, D, LT, p);
}

DECL M casadi_chol(const M& A) {
  return chol(A);
}

DECL M casadi_norm_inf_mul(const M& x, const M& y) {
  return norm_inf_mul(x, y);
}

DECL M casadi_sparsify(const M& A, double tol=0) {
  return sparsify(A, tol);
}

DECL void casadi_expand(const M& ex, M& OUTPUT1, M& OUTPUT2) {
  expand(ex, OUTPUT1, OUTPUT2);
}

DECL M casadi_pw_const(const M &t, const M& tval, const M& val) {
  return pw_const(t, tval, val);
}

DECL M casadi_pw_lin(const M& t, const M& tval, const M& val) {
  return pw_lin(t, tval, val);
}

DECL M casadi_heaviside(const M& x) {
  return heaviside(x);
}

DECL M casadi_rectangle(const M& x) {
  return rectangle(x);
}

DECL M casadi_triangle(const M& x) {
  return triangle(x);
}

DECL M casadi_ramp(const M& x) {
  return ramp(x);
}

DECL M casadi_gauss_quadrature(const M& f, const M& x,
                               const M& a, const M& b,
                               casadi_int order=5) {
  return gauss_quadrature(f, x, a, b, order);
}

DECL M casadi_gauss_quadrature(const M& f, const M& x,
                               const M& a, const M& b,
                               casadi_int order, const M& w) {
  return gauss_quadrature(f, x, a, b, order, w);
}

DECL M casadi_taylor(const M& ex, const M& x, const M& a=0, casadi_int order=1) {
  return taylor(ex, x, a, order);
}

DECL M casadi_mtaylor(const M& ex, const M& x, const M& a, casadi_int order=1) {
  return mtaylor(ex, x, a, order);
}

DECL M casadi_mtaylor(const M& ex, const M& x, const M& a, casadi_int order,
                      const std::vector<casadi_int>& order_contributions) {
  return mtaylor(ex, x, a, order, order_contributions);
}

DECL M casadi_poly_coeff(const M& ex,
                         const M&x) {
  return poly_coeff(ex, x);
}

DECL M casadi_poly_roots(const M& p) {
  return poly_roots(p);
}

DECL M casadi_eig_symbolic(const M& m) {
  return eig_symbolic(m);
}

#endif
%enddef

%define MATRIX_ALL(DECL, FLAG)
MATRIX_FUN(DECL, (FLAG | IS_DMATRIX), Matrix<double>)
MATRIX_FUN(DECL, (FLAG | IS_SX), Matrix<SXElem>)
%enddef

%define MX_FUN(DECL, FLAG, M)
#if FLAG & IS_MEMBER
DECL M casadi_find(const M& x) {
  return find(x);
}
DECL M casadi_low(const M& v, const M& p, const Dict& options = Dict()) {
  return low(v, p, options);
}
DECL M casadi_inv_node(const M& x) {
  return inv_node(x);
}
#endif // FLAG & IS_MEMBER

#if FLAG & IS_GLOBAL
DECL std::vector< M >
casadi_matrix_expand(const std::vector< M >& e,
                     const std::vector< M > &boundary = std::vector< M >(),
                     const Dict& options = Dict()) {
  return matrix_expand(e, boundary, options);
}

DECL M casadi_matrix_expand(const M& e,
                            const std::vector< M > &boundary = std::vector< M >(),
                            const Dict& options = Dict()) {
  return matrix_expand(e, boundary, options);
}

DECL M casadi_graph_substitute(const M& ex, const std::vector< M >& v,
                         const std::vector< M > &vdef) {
  return graph_substitute(ex, v, vdef);
}

DECL std::vector< M >
casadi_graph_substitute(const std::vector< M > &ex,
                 const std::vector< M > &v,
                 const std::vector< M > &vdef) {
  return graph_substitute(ex, v, vdef);
}
DECL M casadi_bspline(const M& x,
        const DM& coeffs,
        const std::vector< std::vector<double> >& knots,
        const std::vector<casadi_int>& degree,
        casadi_int m,
        const Dict& opts = Dict()) {
  return bspline(x, coeffs, knots, degree, m, opts);
}
DECL M casadi_bspline(const M& x,
        const M& coeffs,
        const std::vector< std::vector<double> >& knots,
        const std::vector<casadi_int>& degree,
        casadi_int m,
        const Dict& opts = Dict()) {
  return bspline(x, coeffs, knots, degree, m, opts);
}
DECL MX casadi_bspline(const MX& x,
        const MX& coeffs,
        const std::vector<MX>& knots,
        const std::vector<casadi_int>& degree,
        casadi_int m,
        const Dict& opts = Dict()) {
  return bspline(x, coeffs, knots, degree, m, opts);
}
DECL M casadi_convexify(const M& H,
        const Dict& opts = Dict()) {
  return convexify(H, opts);
}
DECL M casadi_stop_diff(const M& expr, casadi_int order) {
  return stop_diff(expr, order);
}
DECL M casadi_stop_diff(const M& expr, const M& var, casadi_int order) {
  return stop_diff(expr, var, order);
}
DECL std::vector< M > casadi_difference(const std::vector< M >& a, const std::vector< M >& b) {
  return difference(a, b);
}
DECL M casadi_no_hess(const M& expr) {
  return no_hess(expr);
}
DECL M casadi_no_grad(const M& expr) {
  return no_grad(expr);
}
DECL M casadi_kron_contract(const M& m, const M& x, bool inner) {
  return kron_contract(m, x, inner);
}

#endif
%enddef

%define MX_ALL(DECL, FLAG)
MX_FUN(DECL, (FLAG | IS_MX), MX)
%enddef
%include <casadi/core/matrix_fwd.hpp>
%include <casadi/core/matrix_decl.hpp>
%include <casadi/core/dm_fwd.hpp>
%include <casadi/core/sx_fwd.hpp>

// Remove from API
%warnfilter(401) casadi::Matrix<casadi_int>;
%template() casadi::Matrix<casadi_int>;


%template(DM) casadi::Matrix<double>;
%extend casadi::Matrix<double> {
   %template(DM) Matrix<SXElem>;
};


namespace casadi{
  %extend Matrix<double> {
#ifndef SWIGWASMJS
    void assign(const casadi::Matrix<double>&rhs) { (*$self)=rhs; }
#endif
    %matrix_helpers(casadi::Matrix<double>, _MIndex)

  }

}



// Extend DM with SWIG unique features
#if !defined(SWIGWASMJS) && !defined(SWIGJAVA) && !defined(SWIGRUST) && !defined(SWIGDOTNET)
namespace casadi{
  %extend Matrix<double> {
    // Convert to a dense matrix
    GUESTOBJECT* full() const {
      return full(*$self);
    }

    // Convert to a sparse matrix
    GUESTOBJECT* sparse() const {
      return sparse(*$self);
    }
  }

} // namespace casadi
#elif defined(SWIGJAVA) || defined(SWIGRUST)
namespace casadi{
  %extend Matrix<double> {
    // Dense entries, one array per row
    std::vector< std::vector<double> > full() const {
      std::vector<double> d($self->numel());
      casadi::casadi_densify($self->ptr(), $self->sparsity(), casadi::get_ptr(d), false);
      std::vector< std::vector<double> > r($self->size1(), std::vector<double>($self->size2()));
      for (casadi_int i=0; i<$self->size1(); ++i)
        for (casadi_int j=0; j<$self->size2(); ++j) r[i][j] = d[i + j*$self->size1()];
      return r;
    }
  }
} // namespace casadi
#endif


%include <casadi/core/sx_elem.hpp>


namespace casadi {
%extend Matrix<SXElem>{
    %matrix_helpers(casadi::Matrix<casadi::SXElem>, _MIndex)


};

} // namespace casadi


%template(SX) casadi::Matrix<casadi::SXElem>;
%extend casadi::Matrix<casadi::SXElem> {
   %template(SX) Matrix<double>;
};

%include <casadi/core/mx.hpp>

%extend casadi::MX{
  %matrix_helpers(casadi::MX, %arg(_MIndex | MX))
};


#if !defined(SWIGPYTHON) && !defined(SWIGJAVA)
%ignore FunctionBuffer;
%ignore _function_buffer_eval;
#endif

%include <casadi/core/function.hpp>

%include <casadi/core/external.hpp>
%include <casadi/core/integrator.hpp>
%include <casadi/core/conic.hpp>
%include <casadi/core/nlpsol.hpp>
%include <casadi/core/onnx_function.hpp>
%include <casadi/core/graph_model.hpp>
%include <casadi/core/graph_builder.hpp>
%include <casadi/core/rootfinder.hpp>
%include <casadi/core/linsol.hpp>
%include <casadi/core/dple.hpp>
%include <casadi/core/expm.hpp>
%include <casadi/core/interpolant.hpp>
%include <casadi/core/blazing_spline.hpp>
%include <casadi/core/modelica_parser.hpp>

%feature("copyctor", "0") casadi::CodeGenerator;
%include <casadi/core/code_generator.hpp>

#ifdef CASADI_TARGET_MEMBER_FUNS
// Member functions wrapped by the target
CASADI_TARGET_MEMBER_FUNS
#define FLAG IS_GLOBAL
#else
#define FLAG (IS_GLOBAL | IS_MEMBER)
#endif

// Wrap non-member functions, possibly with casadi_ prefix

%inline {
  namespace casadi {
    SPARSITY_INTERFACE_ALL(inline, FLAG)
    GENERIC_EXPRESSION_ALL(inline, FLAG)
    GENERIC_MATRIX_ALL(inline, FLAG)
    MATRIX_ALL(inline, FLAG)
    MX_ALL(inline, FLAG)
  }
}

%inline {
  namespace casadi {
    inline DM casadi_bspline_dual(const std::vector<double>& x,
        const std::vector< std::vector<double> >& knots,
        const std::vector<casadi_int>& degree, const Dict& opts=Dict()) {
      return MX::bspline_dual(x, knots, degree, opts);
    }
  }
}

// Wrap the casadi_ prefixed functions in member functions

%feature("director") casadi::Callback;

// Callback: hand-written stubs for the director-reshaped methods

%include <casadi/core/importer.hpp>
%include <casadi/core/callback.hpp>

%include <casadi/core/global_options.hpp>

%include <casadi/core/casadi_meta.hpp>

%include <casadi/core/integration_tools.hpp>
%include <casadi/core/nlp_tools.hpp>
%include <casadi/core/tools.hpp>
%include <casadi/core/nlp_builder.hpp>
%include <casadi/core/dae_builder.hpp>
%include <casadi/core/xml_file.hpp>
%include <casadi/core/archiver.hpp>
%include <casadi/core/filesystem.hpp>
%include <casadi/core/blas.hpp>
%include <casadi/core/options.hpp>

%feature("copyctor", "0") casadi::SerializerBase;
%feature("copyctor", "0") casadi::DeserializerBase;
%feature("copyctor", "0") casadi::StringSerializer;
%feature("copyctor", "0") casadi::StringDeserializer;
%feature("copyctor", "0") casadi::FileSerializer;
%feature("copyctor", "0") casadi::FileDeserializer;
%nodefaultctor casadi::SerializerBase;
%nodefaultctor casadi::DeserializerBase;


%include <casadi/core/serializer.hpp>


%feature("director") casadi::OptiCallback;


%apply casadi_int &OUTPUT { Opti::ConstraintType &OUTPUT };

%typemap(argout, noblock=1,fragment="casadi_all") casadi::Opti::ConstraintType &OUTPUT {
  %append_output(swig::from(static_cast<casadi_int>(*$1)));
}

%typemap(in, doc="Opti.ConstraintType", pystub_in="Any", noblock=1, numinputs=0) casadi::Opti::ConstraintType &OUTPUT (casadi::Opti::ConstraintType m) {
 $1 = &m;
}


// Zero-argument getters read as properties where the target defines make_property
#ifdef make_property
make_property(casadi::Opti, debug, OptiAdvanced);
make_property(casadi::Opti, advanced, OptiAdvanced);
make_property(casadi::OptiSol, opti, Opti);

%define make_property_opti(name, type)
  make_property(casadi::Opti, name, type);
%enddef

make_property_opti(f, MX)
make_property_opti(g, MX)
make_property_opti(x, MX)
make_property_opti(p, MX)
make_property_opti(lam_g, MX)
make_property_opti(lbg, MX)
make_property_opti(ubg, MX)
make_property_opti(nx, int)
make_property_opti(np, int)
make_property_opti(ng, int)
make_property_opti(x_linear_scale, DM)
make_property_opti(x_linear_scale_offset, DM)
make_property_opti(g_linear_scale, DM)
make_property_opti(f_linear_scale, float)
#endif

%include <casadi/core/optistack.hpp>


%include <casadi/core/resource.hpp>

#ifdef SWIGPYTHON
%include "python/casadi_python_post.i"
#elif defined(SWIGMATLAB)
%include "matlab/casadi_matlab_post.i"
#elif defined(SWIGWASMJS)
%include "wasm-js/casadi_wasmjs_post.i"
#elif defined(SWIGJAVA)
%include "java/casadi_java_post.i"
#endif
// Cleanup for dependent modules
%exception {
  $action
}
