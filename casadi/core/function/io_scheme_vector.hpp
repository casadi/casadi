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
    // Data members (all public)
  public:
    /// Vector of data
    std::vector<T> data;

    /// Scheme
    std::vector<std::string> str;

  public:
    IOSchemeVector(const std::vector<T>& d = std::vector<T>(),
                   const std::vector<std::string>& s = std::vector<std::string>()) : data(d), str(s) {}
    IOSchemeVector(const std::vector<T>& d, InputOutputScheme scheme) : data(d), str(IOScheme(scheme)) {}

#ifndef SWIG
    T operator[](int i) const {
      casadi_assert_message(i>=0 && i<this->data.size(), "Index error: supplied integer must be >=0 and <= " << this->data.size() << " but got " << i << ".");
      return this->data.at(i);
    }
    T operator[](const std::string& name) const {
      for (std::vector<std::string>::const_iterator i=str.begin(); i!=str.end(); ++i) {
        size_t col = i->find(':');
        if (i->compare(0, col, name)==0) {
          return (*this)[i-str.begin()];
        }
      }
      std::stringstream entries;
      for (size_t i=0; i!=str.size(); ++i) {
         if (i!=0) entries << ", ";
         size_t col = str[i].find(':');
         entries << str[i].substr(0, col);
      }
      casadi_error("entry '" << name << "' not available. Available entries are " << entries.str());
      return T();
    }
#endif // SWIG
#ifndef SWIGMATLAB
    T __getitem__(int i) const { if (i<0) i+= this->data.size(); return (*this)[i]; }
    T __getitem__(const std::string& name) const { return (*this)[name]; }
    int __len__() const { return this->data.size(); }
#endif // SWIGMATLAB
    std::vector<T> vector() const { return this->data; }

    /// Print a description of the object
    void print(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const {
      stream << "IOSchemeVector(" ;
      for (int i=0; i<this->data.size(); ++i) {
        if (this->str.empty()) {
          stream << "none";
        } else {
          str.at(i);
        }
        stream << "=" << this->data[i];
        if (i<this->data.size()-1) stream << ", ";
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
