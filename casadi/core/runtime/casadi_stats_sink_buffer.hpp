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

// The sink of a stats stream without function pointers (codegen option
// allow_function_pointers off): a plain counter, bumped under the owner's mutex if it gave one.
// mutex is a CASADI_MUTEX_TYPE*, of the CASADI_THREAD_TYPE everyone writing to the sink is
// compiled with; declared void*, so that this file needs no threads runtime.
// Declarations only: generated code emits this file into its header as well.

// FILTER-MACROS OFF
#ifndef CASADI_STATS_ROOT
#define CASADI_STATS_ROOT (-1)
#endif
#ifndef CASADI_STATS_STRUCT
#define CASADI_STATS_STRUCT 2
struct casadi_stats_sink {
  unsigned char* p;
  casadi_int cap;
  casadi_int needed;
  void* mutex;
};
#elif CASADI_STATS_STRUCT != 2
#error "struct casadi_stats_sink: code generated with and without allow_function_pointers"
#endif
// FILTER-MACROS ON
