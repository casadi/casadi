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
#include <map>
namespace casadi {

/** \brief A vector container with associated IOScheme
  A class
*/
  template<typename T>
  class CASADI_EXPORT IOSchemeVector : public PrintableObject<IOSchemeVector<T> > {
  public:
    std::pair<std::map<std::string, T>, std::vector<std::string> > v3;

    IOSchemeVector(const std::map<std::string, T>& m,
                   const std::vector<std::string>& s) : v3(make_pair(m, s)) {
    }

    IOSchemeVector(const std::vector<T>& d = std::vector<T>(),
                   const std::vector<std::string>& s = std::vector<std::string>()) {
      // Get/form string
      if (s.empty()) {
        v3.second.resize(d.size());
        for (size_t i=0; i!=d.size(); ++i) {
          std::stringstream ss;
          ss << i;
          v3.second[i] = ss.str();
        }
      } else {
        casadi_assert(d.size()==s.size());
        v3.second = s;
      }

      // Create map
      for (size_t i=0; i!=d.size(); ++i) {
        v3.first[v3.second[i]] = d[i];
      }
    }

    IOSchemeVector(const std::vector<T>& d, InputOutputScheme scheme) {
      v3.second = IOScheme(scheme);
      casadi_assert(d.size()==v3.second.size());
      for (size_t i=0; i!=d.size(); ++i) {
        v3.first[v3.second[i]] = d[i];
      }
    }

#ifndef SWIG
    T operator[](int i) const {
      casadi_assert_message(i>=0 && i<this->v3.second.size(), "Index error: supplied integer must be >=0 and <= " << this->v3.second.size() << " but got " << i << ".");
      std::string s = this->v3.second.at(i);
      typename std::map<std::string, T>::const_iterator it;
      it=this->v3.first.find(s);
      casadi_assert(it!=this->v3.first.end());
      return it->second;
    }
#endif // SWIG
#ifndef SWIGMATLAB
    T __getitem__(int i) const { if (i<0) i+= this->v3.second.size(); return (*this)[i]; }
    int __len__() const { return this->v3.second.size(); }
#endif // SWIGMATLAB

    /// Print a description of the object
    void print(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const {
      stream << "IOSchemeVector(" ;
      for (int i=0; i<this->v3.second.size(); ++i) {
        if (this->v3.second.at(i).empty()) {
          stream << "none";
        } else {
          stream << this->v3.second.at(i);
        }
        typename std::map<std::string, T>::const_iterator it;
        it=this->v3.first.find(this->v3.second.at(i));
        casadi_assert(it!=this->v3.first.end());
        stream << "=" << it->second;
        if (i<this->v3.second.size()-1) stream << ", ";
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
