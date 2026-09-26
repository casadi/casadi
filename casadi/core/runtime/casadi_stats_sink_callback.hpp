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

// Stats sink: buffer p of cap bytes, needed bytes reserved so far; optional owner reserve,
// and commit, called once the record reserved at pos is written

// FILTER-MACROS OFF
#ifndef CASADI_STATS_SINK
#define CASADI_STATS_SINK
struct casadi_stats_sink {
  unsigned char* p;
  casadi_int cap;
  casadi_int needed;
  unsigned char* (*reserve)(struct casadi_stats_sink* s, casadi_int len, casadi_int* pos);
  void (*commit)(struct casadi_stats_sink* s, casadi_int pos);
  void* data;
};
#endif
// FILTER-MACROS ON
