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


#ifndef CASADI_PROFILING_HPP
#define CASADI_PROFILING_HPP

#include <fstream>
#include <cstring>
#include <iostream>

#include "casadi_common.hpp"

namespace casadi {
/// \cond INTERNAL

/*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */

/**
 * Returns the real time, in seconds, or -1.0 if an error occurred.
 *
 * Time is measured since an arbitrary and OS-dependent start time.
 * The returned real time is only useful for computing an elapsed time
 * between two calls to this function.
 */
CASADI_CORE_EXPORT double getRealTime();

enum ProfilingData_Type { ProfilingData_Type_TIMELINE,
                          ProfilingData_Type_SOURCE,
                          ProfilingData_Type_NAME,
                          ProfilingData_Type_ENTRY,
                          ProfilingData_Type_EXIT,
                          ProfilingData_Type_IO };

enum ProfilingData_FunctionType { ProfilingData_FunctionType_MXFunction,
                                  ProfilingData_FunctionType_SXFunction,
                                  ProfilingData_FunctionType_Other };

struct ProfilingHeader {
  ProfilingData_Type type;
};

struct ProfilingData_TIMELINE {
  double local;
  double total;
  long thisp;
  int line_number;
};

struct ProfilingData_SOURCE {
  long thisp;
  int line_number;
  int length;
  int opcode;
  long dependency;
};

struct ProfilingData_NAME {
  long thisp;
  int length;
  ProfilingData_FunctionType type;
  int algorithm_size;
  int numin;
  int numout;
};

struct ProfilingData_IO {
  int nrow;
  int ncol;
  int ndata;
};

struct ProfilingData_ENTRY {
  long thisp;
};

struct ProfilingData_EXIT {
  double total;
  long thisp;
};

template<typename T>
ProfilingData_Type ProfilingType();

template<>
inline ProfilingData_Type ProfilingType<ProfilingData_TIMELINE>()
{ return ProfilingData_Type_TIMELINE; }
template<>
inline ProfilingData_Type ProfilingType<ProfilingData_SOURCE>()
{ return ProfilingData_Type_SOURCE; }
template<>
inline ProfilingData_Type ProfilingType<ProfilingData_NAME>()
{ return ProfilingData_Type_NAME; }
template<>
inline ProfilingData_Type ProfilingType<ProfilingData_ENTRY>()
{ return ProfilingData_Type_ENTRY; }
template<>
inline ProfilingData_Type ProfilingType<ProfilingData_EXIT>()
{ return ProfilingData_Type_EXIT; }
template<>
inline ProfilingData_Type ProfilingType<ProfilingData_IO>()
{ return ProfilingData_Type_IO; }


template<typename T>
void profileWrite(std::ofstream &f, const T& s) {
  ProfilingHeader hd;
  hd.type   = ProfilingType<T>();
  f.write(reinterpret_cast<const char*>(&hd), sizeof(hd));
  f.write(reinterpret_cast<const char*>(&s), sizeof(s));
}

template<typename T>
void profileWriteBare(std::ofstream &f, const T& s) {
  f.write(reinterpret_cast<const char*>(&s), sizeof(s));
}

template<typename T>
long ptrToLong(T *a) {
  long r;
  std::memcpy(&r, &a, sizeof(T*));
  return r;
}

template<typename T>
void profileWriteName(std::ofstream &f, T *a, const std::string &name,
                      ProfilingData_FunctionType type, int algorithm_size) {
  ProfilingData_NAME s;
  s.thisp=ptrToLong(a);
  s.length=name.size();
  s.type = type;
  s.algorithm_size = algorithm_size;
  s.numin = a->getNumInputs();
  s.numout = a->getNumOutputs();
  profileWrite(f, s);
  f << name;
  for (int i=0;i<s.numin;++i) {
    ProfilingData_IO ss;
    ss.nrow = a->input(i).size1();
    ss.ncol = a->input(i).size2();
    ss.ndata = a->input(i).size();
    profileWriteBare(f, ss);
  }
  for (int i=0;i<s.numout;++i) {
    ProfilingData_IO ss;
    ss.nrow = a->output(i).size1();
    ss.ncol = a->output(i).size2();
    ss.ndata = a->output(i).size();
    profileWriteBare(f, ss);
  }
}

template<typename T>
void profileWriteEntry(std::ofstream &f, T *a) {
  ProfilingData_ENTRY s;
  s.thisp=ptrToLong(a);
  profileWrite(f, s);
}

template<typename T>
void profileWriteExit(std::ofstream &f, T *a, double total) {
  ProfilingData_EXIT s;
  s.thisp=ptrToLong(a);
  s.total=total;
  profileWrite(f, s);
}

template<typename T>
void profileWriteTime(std::ofstream &f, T *a, int line_number, double local, double total) {
  ProfilingData_TIMELINE s;
  s.local = local;
  s.total = total;
  s.thisp = ptrToLong(a);
  s.line_number = line_number;
  profileWrite(f, s);
}


template<typename T, typename T2>
void profileWriteSourceLine(std::ofstream &f, T *a, int line_number, const std::string &sourceline,
                            int opcode, T2 *dependency) {
  ProfilingData_SOURCE s;
  s.thisp = ptrToLong(a);
  s.line_number = line_number;
  s.length = sourceline.size();
  s.opcode = opcode;
  s.dependency = ptrToLong(dependency);
  profileWrite(f, s);
  f << sourceline;
}

template<typename T>
void profileWriteSourceLine(std::ofstream &f, T *a, int line_number,
                            const std::string &sourceline, int opcode) {
  ProfilingData_SOURCE s;
  s.thisp = ptrToLong(a);
  s.line_number = line_number;
  s.length = sourceline.size();
  s.opcode = opcode;
  s.dependency = 0;
  profileWrite(f, s);
  f << sourceline;
}

/// \endcond
} // namespace casadi

#endif // CASADI_PROFILING_HPP
