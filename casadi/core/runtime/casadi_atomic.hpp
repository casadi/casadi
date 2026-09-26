//
//    MIT No Attribution
//
//    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//    Permission is hereby granted, free of charge, to any person obtaining a copy of this
//    software and associated documentation files (the "Software"), to deal in the Software
//    without restriction, including without limitation the rights to use, copy, modify,
//    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//    permit persons to whom the Software is furnished to do so.
//
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//

// FILTER-MACROS OFF
/* Relaxed atomic fetch-and-add on an integer; undefined if unsupported (mutex fallback) */
#ifndef CASADI_ATOMIC_FETCH_ADD
  #if defined(CASADI_THREAD_TYPE) && CASADI_THREAD_TYPE == CASADI_THREAD_TYPE_NONE
    #define CASADI_ATOMIC_FETCH_ADD(p, n) ((*(p) += (n)) - (n))
  #elif defined(_MSC_VER) && !defined(__clang__) && (defined(_M_X64) || defined(_M_ARM64))
    #include <intrin.h>
    #define CASADI_ATOMIC_FETCH_ADD(p, n) (sizeof(*(p)) == 8 ? \
      _InterlockedExchangeAdd64((volatile long long*) (p), (long long) (n)) : \
      _InterlockedExchangeAdd((volatile long*) (p), (long) (n)))
  #elif defined(__GNUC__) || defined(__clang__)
    #define CASADI_ATOMIC_FETCH_ADD(p, n) __atomic_fetch_add(p, n, __ATOMIC_RELAXED)
  #endif
#endif

// FILTER-MACROS ON
