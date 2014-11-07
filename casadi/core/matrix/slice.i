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

%include <casadi/core/printable_object.i>

%include <casadi/core/matrix/slice.hpp>

namespace casadi{

%{
#ifndef SWIGXML

/// casadi::Slice
template<> char meta< casadi::Slice >::expected_message[] = "Expecting Slice or number";
template <>
int meta< casadi::Slice >::as(GUESTOBJECT *p, casadi::Slice &m) {
  NATIVERETURN(casadi::Slice,m)
#ifdef SWIGPYTHON
  if (PyInt_Check(p)) {
    m.start_ = PyInt_AsLong(p);
    m.stop_ = m.start_+1;
    if (m.stop_==0) m.stop_ = std::numeric_limits<int>::max();
    return true;
  } else if (PySlice_Check(p)) {
    PySliceObject *r = (PySliceObject*)(p);
    m.start_ = (r->start == Py_None || PyInt_AsLong(r->start) < std::numeric_limits<int>::min()) ? std::numeric_limits<int>::min() : PyInt_AsLong(r->start);
    m.stop_  = (r->stop ==Py_None || PyInt_AsLong(r->stop)> std::numeric_limits<int>::max()) ? std::numeric_limits<int>::max() : PyInt_AsLong(r->stop) ;
    if(r->step !=Py_None) m.step_  = PyInt_AsLong(r->step);
    return true;
  } else {
    return false;
  }
#else
  return false;
#endif
}

/// casadi::IndexList
template<> char meta< casadi::IndexList >::expected_message[] = "Expecting Slice or number or list of ints";
template <>
int meta< casadi::IndexList >::as(GUESTOBJECT *p, casadi::IndexList &m) {
  NATIVERETURN(casadi::IndexList,m)
  if (meta< int >::couldbe(p)) {
    m.type = casadi::IndexList::INT;
    meta< int >::as(p,m.i);
  } else if (meta< std::vector<int> >::couldbe(p)) {
    m.type = casadi::IndexList::IVECTOR;
    return meta< std::vector<int> >::as(p,m.iv);
  } else if (meta< casadi::Slice>::couldbe(p)) {
    m.type = casadi::IndexList::SLICE;
    return meta< casadi::Slice >::as(p,m.slice);
  } else {
    return false;
  }
  return true;
}

#endif // SWIGXML

%}

%{
template<> swig_type_info** meta< casadi::Slice >::name = &SWIGTYPE_p_casadi__Slice;
template<> swig_type_info** meta< casadi::IndexList >::name = &SWIGTYPE_p_casadi__IndexList;
%}

%my_generic_const_typemap(PRECEDENCE_SLICE,casadi::Slice);
%my_generic_const_typemap(PRECEDENCE_IndexVector,casadi::IndexList);

} // namespace casadi

#endif // CASADI_SLICE_I
