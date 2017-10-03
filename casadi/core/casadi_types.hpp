/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
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


#ifndef CASADI_CASADI_TYPES_HPP
#define CASADI_CASADI_TYPES_HPP

#include <climits>
#include <vector>
#include <utility>

#include "casadi_logger.hpp"

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

// warinng C4251: Need a dll interface?
#pragma warning(disable:4251)

// warning C4715: Not all control paths return a value
#pragma warning(disable:4715)

// warning C4800: 'int' : forcing value to bool 'true'or 'false'(performance warning)
#pragma warning(disable:4800)

// warning C4910: __declspec(dllexport) and extern incompatible on an explicit instantiation
#pragma warning(disable:4910)

// ?
#pragma warning(disable:4996)

#endif // _MSC_VER

  // Macro "minor" is sometimes defined, cf.
  // https://stackoverflow.com/questions/22240973/major-and-minor-macros-defined-in-sys-sysmacros-h-pulled-in-by-iterator
#undef minor

  // The number of derivative directions for which the tool has been optimized
  const int optimized_num_dir = 64;

  // Type with a size corresponding to that of double (or smaller) that can be used to hold a set
  // of booleans. If the compiler supports C99 or has defined __SIZEOF_LONG_LONG__,
  // we shall use the long long datatype, which is 64 bits, otherwise long
  #if (defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L || defined(__SIZEOF_LONG_LONG__))
  typedef unsigned long long bvec_t;
  #else
  typedef unsigned long bvec_t;
  #endif

  // Number of directions we can deal with at a time
  // the size of bvec_t in bits (CHAR_BIT is the number of bits per byte, usually 8)
  const int bvec_size = CHAR_BIT*sizeof(bvec_t);

  // Make sure that the integer datatype is indeed smaller or equal to the double
  //assert(sizeof(bvec_t) <= sizeof(double)); // doesn't work - very strange

  ///@{
  /** \brief  Function pointer types for the C API */
  typedef void (*signal_t)(void);
  typedef int (*getint_t)(void);
  typedef const char* (*name_t)(int i);
  typedef const int* (*sparsity_t)(int i);
  typedef void* (*alloc_mem_t)(void);
  typedef int (*init_mem_t)(void* mem);
  typedef void (*free_mem_t)(void* mem);
  typedef int (*work_t)(int* sz_arg, int* sz_res, int* sz_iw, int* sz_w);
  typedef int (*eval_t)(const double** arg, double** res, int* iw, double* w, void* mem);
  ///@}

#endif // SWIG

} // namespace casadi

#endif // CASADI_CASADI_TYPES_HPP
