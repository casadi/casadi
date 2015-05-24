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


#include "io_scheme_internal.hpp"
#include <map>
using namespace std;

namespace casadi {

    IOSchemeCustomInternal::IOSchemeCustomInternal(const std::vector<std::string>& data) :
      data_(data) {
    }

    std::string IOSchemeCustomInternal::name() const {
      return "customIO";
    }

    std::string IOSchemeCustomInternal::entryNames() const {
      std::stringstream ss;
      for (int i=0; i<data_.size();++i) {
         if (i!=0) ss << ", ";
         size_t col = data_[i].find(':');
         ss << data_[i].substr(0, col);
      }
      return ss.str();
    }

    std::string IOSchemeCustomInternal::entry(int i) const {
      casadi_assert_message(i>=0 && i<data_.size(),
                            "customIO::entry(): requesting entry for index " << i
                            << ", but IOScheme is only length " << data_.size());
      size_t col = data_[i].find(':');
      return data_[i].substr(0, col);
    }

    std::string IOSchemeCustomInternal::entryEnum(int i) const {
      return "";
    }

    std::string IOSchemeCustomInternal::describe(int i) const {
      casadi_assert_message(i>=0 && i<data_.size(),
                            "customIO::entry(): requesting entry for index " << i
                            << ", but IOScheme is only length " << data_.size());
      return data_[i];
    }

    int IOSchemeCustomInternal::index(const std::string &name) const {
      for (vector<string>::const_iterator i=data_.begin(); i!=data_.end(); ++i) {
        size_t col = i->find(':');
        if (i->compare(0, col, name)==0) return i-data_.begin();
      }
      casadi_error("customIO::index(): entry '" << name
                   << "' not available. Available entries are " << entryNames());
      return -1;
    }

    int IOSchemeCustomInternal::size() const {
      return data_.size();
    }

    void IOSchemeCustomInternal::print(std::ostream &stream) const {
      stream << "customIO(" << entryNames() << ")";
    }

    void IOSchemeCustomInternal::repr(std::ostream &stream) const {
      stream << "customIO(" << entryNames() << ")";
    }

} // namespace casadi
