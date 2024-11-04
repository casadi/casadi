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


#ifndef CASADI_COMMON_HPP
#define CASADI_COMMON_HPP

#include <algorithm>
#include <climits>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#ifdef SWIG
#define SWIG_IF_ELSE(is_swig, not_swig) is_swig
#define SWIG_OUTPUT(arg) OUTPUT
#define SWIG_INOUT(arg) INOUT
#define SWIG_CONSTREF(arg) const arg
#ifdef SWIGMATLAB
#define SWIG_IND1 true
#else // SWIGMATLAB
#define SWIG_IND1 false
#endif // SWIGMATLAB
#else // SWIG
#define SWIG_IF_ELSE(is_swig, not_swig) not_swig
#define SWIG_OUTPUT(arg) arg
#define SWIG_INOUT(arg) arg
#define SWIG_CONSTREF(arg) const arg &
#define SWIG_IND1 false
#endif // SWIG

#include "casadi_types.hpp"

namespace casadi {

  /// Forward declarations
  class SXElem;
  class MX;
  template<class T> class Matrix;
  class Function;
  class Sparsity;
  class CodeGenerator;
  class NlpBuilder;
  struct Variable;
  class DaeBuilder;
  class XmlFile;
  class Importer;

#ifndef SWIG
// Get GCC version if GCC is used
#ifdef __GNUC__
#ifdef __GNUC_MINOR__
#ifdef __GNUC_PATCHLEVEL__
#define GCC_VERSION (__GNUC__ * 10000 +__GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#endif // __GNUC_PATCHLEVEL__
#endif // __GNUC_MINOR__
#endif // __GNUC__

// Disable some Visual studio warnings
#ifdef _MSC_VER

// warning C4018: '<' : signed/unsigned mismatch
#pragma warning(disable:4018)

// warning C4244: Potential loss of data converting double to int
#pragma warning(disable:4244)

// warning C4251: Need a dll interface?
#pragma warning(disable:4251)

// warning C4275: non dll-interface class 'std::exception' used as base for dll-interface
// class 'casadi::CasadiException'
#pragma warning(disable:4275)

// warning C4715: Not all control paths return a value
#pragma warning(disable:4715)

// warning C4800: 'int' : forcing value to bool 'true'or 'false'(performance warning)
#pragma warning(disable:4800)

// warning C4910: __declspec(dllexport) and extern incompatible on an explicit instantiation
#pragma warning(disable:4910)

// warning C4996: 'sprintf': This function or variable may be unsafe. Consider using sprintf_s
// instead
#pragma warning(disable:4996)

#endif // _MSC_VER

  // Macro "minor" is sometimes defined, cf.
  // https://stackoverflow.com/questions/22240973/major-and-minor-macros-defined-in-sys-sysmacros-h-pulled-in-by-iterator
#undef minor

  // Type with a size corresponding to that of double (or smaller) that can be used to hold a set
  // of booleans
  typedef unsigned long long bvec_t;

  // Number of directions we can deal with at a time
  // the size of bvec_t in bits (CHAR_BIT is the number of bits per byte, usually 8)
  const int bvec_size = CHAR_BIT*sizeof(bvec_t);

  // Make sure that the integer datatype is indeed smaller or equal to the double
  //assert(sizeof(bvec_t) <= sizeof(double)); // doesn't work - very strange

  ///@{
  /** \brief  Function pointer types for the C API

      \identifier{7j} */
  typedef int (*config_t)(int, const char**);
  typedef void (*signal_t)(void);
  typedef casadi_int (*getint_t)(void);
  typedef double (*default_t)(casadi_int i);
  typedef const char* (*name_t)(casadi_int i);
  typedef const casadi_int* (*sparsity_t)(casadi_int i);
  typedef int (*diff_t)(casadi_int i);
  typedef int (*casadi_checkout_t)(void);
  typedef void (*casadi_release_t)(int);
  typedef int (*work_t)(casadi_int* sz_arg, casadi_int* sz_res,
    casadi_int* sz_iw, casadi_int* sz_w);
  typedef int (*eval_t)(const double** arg, double** res,
                        casadi_int* iw, double* w, int);
  ///@}

  // Easier to maintain than an enum (serialization/codegen)
  constexpr casadi_int LOOKUP_LINEAR = 0;
  constexpr casadi_int LOOKUP_EXACT = 1;
  constexpr casadi_int LOOKUP_BINARY = 2;

  /// String representation, any type
  template<typename T>
  std::string str(const T& v);

  /// String representation, CasADi type
  template<typename T>
  std::string str(const T& v, bool more);

  /// String representation of vector
  template<typename T>
  std::string str(const std::vector<T>& v, bool more=false);

  /// String representation of set
  template<typename T>
  std::string str(const std::set<T>& v, bool more=false);

  /// String representation of pair
  template<typename T1, typename T2>
  std::string str(const std::pair<T1, T2>& p, bool more=false);

  /// String representation of a map
  template<typename T1, typename T2>
  std::string str(const std::map<T1, T2> &p, bool more=false);

  /// String representation of a dictionary
  template<typename T2>
  std::string str(const std::map<std::string, T2> &p, bool more=false);

  /// String representation of an array
  template<typename T, size_t N>
  std::string str(const std::array<T, N> &p, bool more=false);

  //! \brief Create a list of strings from __VA_ARGS__, no argument
  inline std::vector<std::string> strvec() {
    return {};
  }

  //! \brief Create a list of strings from __VA_ARGS__, one argument
  template<typename T1>
  inline std::vector<std::string> strvec(const T1& t1) {
    return {str(t1)};
  }

  //! \brief Create a list of strings from __VA_ARGS__, two arguments
  template<typename T1, typename T2>
  inline std::vector<std::string> strvec(const T1& t1, const T2& t2) {
    return {str(t1), str(t2)};
  }

  //! \brief Create a list of strings from __VA_ARGS__, three arguments
  template<typename T1, typename T2, typename T3>
  inline std::vector<std::string> strvec(const T1& t1, const T2& t2, const T3& t3) {
    return {str(t1), str(t2), str(t3)};
  }

  //! \brief Create a list of strings from __VA_ARGS__, four arguments
  template<typename T1, typename T2, typename T3, typename T4>
  inline std::vector<std::string> strvec(const T1& t1, const T2& t2, const T3& t3,
                                         const T4& t4) {
    return {str(t1), str(t2), str(t3), str(t4)};
  }

  //! \brief Create a list of strings from __VA_ARGS__, five arguments
  template<typename T1, typename T2, typename T3, typename T4, typename T5>
  inline std::vector<std::string> strvec(const T1& t1, const T2& t2, const T3& t3,
                                         const T4& t4, const T5& t5) {
    return {str(t1), str(t2), str(t3), str(t4), str(t5)};
  }

  //! \brief Create a list of strings from __VA_ARGS__, six arguments
  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  inline std::vector<std::string> strvec(const T1& t1, const T2& t2, const T3& t3,
                                         const T4& t4, const T5& t5, const T6& t6) {
    return {str(t1), str(t2), str(t3), str(t4), str(t5), str(t6)};
  }

  //! \brief Create a string from a formatted string
  inline std::string fmtstr(const std::string& fmt, const std::vector<std::string>& args) {
    std::string s = fmt;
    for (auto&& e : args) {
      std::string::size_type n = s.find("%s");
      if (n==std::string::npos) return "** Ill-formatted string ** " + fmt;
      s.replace(n, 2, e);
    }
    return s;
  }

  // Implementations
  template<typename T>
  std::string str(const T& v) {
    std::stringstream ss;
    ss << v;
    return ss.str();
  }

  template<typename T>
  std::string str(const T& v, bool more) {
    return v.get_str(more);
  }

  template<typename T>
  std::string str(const std::vector<T>& v, bool more) {
    std::stringstream ss;
    ss << "[";
    for (casadi_int i=0; i<v.size(); ++i) {
      if (i!=0) ss << ", ";
      ss << v[i];
    }
    ss << "]";
    return ss.str();
  }

  template<typename T>
  std::string str(const std::set<T>& v, bool more) {
    std::stringstream ss;
    ss << "{";
    casadi_int cnt = 0;
    for (const auto& e : v) {
      if (cnt++!=0) ss << ", ";
      ss << e;
    }
    ss << "}";
    return ss.str();
  }

  template<typename T1, typename T2>
  std::string str(const std::pair<T1, T2>& p, bool more) {
    std::stringstream ss;
    ss << "[" << p.first << "," << p.second << "]";
    return ss.str();
  }

  template<typename T1, typename T2>
  std::string str(const std::map<T1, T2>& p, bool more) {
    std::stringstream ss;
    ss << "{";
    casadi_int count = 0;
    for (auto& e : p) {
      ss << e.first << ": " << e.second;
      if (++count < p.size()) ss << ", ";
    }
    ss << "}";
    return ss.str();
  }

  template<typename T2>
  std::string str(const std::map<std::string, T2>& p, bool more) {
    std::stringstream ss;
    ss << "{";
    casadi_int count = 0;
    for (auto& e : p) {
      ss << "\"" << e.first << "\": " << e.second;
      if (++count < p.size()) ss << ", ";
    }
    ss << "}";
    return ss.str();
  }

  template<typename T, size_t N>
  std::string str(const std::array<T, N> &v, bool more) {
    std::stringstream ss;
    ss << "[";
    for (casadi_int i=0; i<N; ++i) {
      if (i!=0) ss << ", ";
      ss << v[i];
    }
    ss << "]";
    return ss.str();
  }

#endif // SWIG

  template<typename _Mutex>
  class conditional_lock_guard {
  public:
    typedef _Mutex mutex_type;

    conditional_lock_guard(mutex_type& m, bool condition) : mtx_(m), condition_(condition) {
      if (condition_) mtx_.lock();
    }

    ~conditional_lock_guard() {
      if (condition_) mtx_.unlock();
    }

    conditional_lock_guard(const conditional_lock_guard&) = delete;
    conditional_lock_guard& operator=(const conditional_lock_guard&) = delete;

  private:
    mutex_type&  mtx_;
    bool condition_;
  };

} // namespace casadi

#include "casadi_logger.hpp"

#endif // CASADI_COMMON_HPP
