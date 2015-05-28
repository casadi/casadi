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


#ifndef CASADI_IO_SCHEME_VECTOR_HPP
#define CASADI_IO_SCHEME_VECTOR_HPP

#include "schemes_metadata.hpp"
#include "io_scheme.hpp"
#include "../casadi_exception.hpp"
#include "../printable_object.hpp"
namespace casadi {

/** \brief A vector container with associated IOScheme
  A class
*/
  template<typename T>
  class CASADI_EXPORT IOSchemeVector : public PrintableObject<IOSchemeVector<T> > {
  public:
    std::pair<std::vector<T>, std::vector<std::string> > v2;

    IOSchemeVector(const std::map<std::string, T>& m,
                   const std::vector<std::string>& s) {
      v2.second = s;
      v2.first.resize(s.size());
      for (size_t i=0; i!=s.size(); ++i) {
        typename std::map<std::string, T>::const_iterator it=m.find(s[i]);
        if (it!=m.end()) {
          v2.first[i]=it->second;
        }
      }
    }

    IOSchemeVector(const std::vector<T>& d = std::vector<T>(),
                   const std::vector<std::string>& s = std::vector<std::string>()) {
      v2.first = d;
      v2.second = s;
      if (v2.second.empty()) {
        v2.second.resize(v2.first.size());
      }
      casadi_assert(v2.first.size()==v2.second.size());
    }
    IOSchemeVector(const std::vector<T>& d, InputOutputScheme scheme) {
      v2.first = d;
      v2.second = IOScheme(scheme);
      casadi_assert(v2.second.size()==v2.first.size());
    }

#ifndef SWIG
    T operator[](int i) const {
      casadi_assert_message(i>=0 && i<this->v2.first.size(), "Index error: supplied integer must be >=0 and <= " << this->v2.first.size() << " but got " << i << ".");
      return this->v2.first.at(i);
    }
#endif // SWIG
#ifndef SWIGMATLAB
    T __getitem__(int i) const { if (i<0) i+= this->v2.first.size(); return (*this)[i]; }
    int __len__() const { return this->v2.first.size(); }
#endif // SWIGMATLAB
    std::vector<T> vector() const {
      return v2.first;
    }

    /// Print a description of the object
    void print(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const {
      stream << "IOSchemeVector(" ;
      for (int i=0; i<this->v2.first.size(); ++i) {
        if (this->v2.second.at(i).empty()) {
          stream << "none";
        } else {
          stream << this->v2.second.at(i);
        }
        stream << "=" << this->v2.first[i];
        if (i<this->v2.first.size()-1) stream << ", ";
      }

      stream <<  ")";
      if (trailing_newline) stream << std::endl;
    }

    /// Print a representation of the object
    void repr(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const {
      print(stream, trailing_newline);
    }
  };


} // namespace casadi


#endif // CASADI_IO_SCHEME_VECTOR_HPP
