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
#ifndef CASADI_SLICE_I
#define CASADI_SLICE_I

%include "printable_object.i"

%include <casadi/core/matrix/slice.hpp>

namespace casadi{

#ifndef SWIGXML

%fragment("to"{Slice}, "header", fragment="fwd") {
   int to_Slice(GUESTOBJECT *p, void *mv, int offs) {
    casadi::Slice *m = static_cast<casadi::Slice*>(mv);
    if (m) m += offs;
    // CasADi-Slice already
    if (is_a(p, $descriptor(casadi::Slice *))) {
      casadi::Slice *mp;
      if (SWIG_ConvertPtr(p, (void **) &mp, $descriptor(casadi::Slice *), 0) == -1)
        return false;
      if(m) *m = *mp;
      return true;
    }
#ifdef SWIGPYTHON
    // Python int
    if (PyInt_Check(p)) {
      if (m) {
        m->start_ = PyInt_AsLong(p);
        m->stop_ = m->start_+1;
        if (m->stop_==0) m->stop_ = std::numeric_limits<int>::max();
      }
      return true;
    }
    // Python slice
    if (PySlice_Check(p)) {
      PySliceObject *r = (PySliceObject*)(p);
      if (m) {
        m->start_ = (r->start == Py_None || PyInt_AsLong(r->start) < std::numeric_limits<int>::min()) 
          ? std::numeric_limits<int>::min() : PyInt_AsLong(r->start);
        m->stop_  = (r->stop ==Py_None || PyInt_AsLong(r->stop)> std::numeric_limits<int>::max())
          ? std::numeric_limits<int>::max() : PyInt_AsLong(r->stop);
        if(r->step !=Py_None) m->step_  = PyInt_AsLong(r->step);
      }
      return true;
    }
#endif // SWIGPYTHON
#ifdef SWIGMATLAB
    if (mxIsChar(p) && mxGetM(p)==1 && mxGetN(p)==1) {
      char ch;
      if(mxGetString(p, &ch,(mwSize)sizeof(&ch))) return SWIG_TypeError;
      if (ch==':') {
        if (m) *m = casadi::Slice();
        return true;
      }
    }
#endif // SWIGMATLAB
    // Failure if reached this point
    return false;
  }
}
%casadi_typemaps_constref(Slice, PRECEDENCE_SLICE, casadi::Slice)

#endif // SWIGXML

} // namespace casadi

#endif // CASADI_SLICE_I
