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

    IOSchemeBuiltinInternal::IOSchemeBuiltinInternal(InputOutputScheme scheme) : scheme_(scheme) {}

    std::string IOSchemeBuiltinInternal::name() const {
      return getSchemeName(scheme_);
    }

    std::string IOSchemeBuiltinInternal::entryNames() const {
      return getSchemeEntryNames(scheme_);
    }

    std::string IOSchemeBuiltinInternal::entry(int i) const {
      return getSchemeEntryName(scheme_, i);
    }

    std::string IOSchemeBuiltinInternal::entryEnum(int i) const {
      return getSchemeEntryEnumName(scheme_, i);
    }

    std::string IOSchemeBuiltinInternal::describe(int i) const {
      std::stringstream ss;
      ss << entryEnum(i) <<  " aka '" << entry(i) << "'";
      return ss.str();
    }

    int IOSchemeBuiltinInternal::index(const std::string &name) const {
      return getSchemeEntryEnum(scheme_, name);
    }

    int IOSchemeBuiltinInternal::size() const {
      return getSchemeSize(scheme_);
    }

    void IOSchemeBuiltinInternal::print(std::ostream &stream) const {
      stream << "builtinIO(" << name() << ")";
    }

    void IOSchemeBuiltinInternal::repr(std::ostream &stream) const {
      stream << "builtinIO(" << name() << ")";
    }

    IOSchemeCustomInternal::IOSchemeCustomInternal(const std::vector<std::string> &entries,
                                                   const std::vector<std::string> &descriptions) :
        entries_(entries), descriptions_(descriptions)  {
      for (int i=0;i<entries.size();++i) {
        entrymap_[entries[i]] = i;
      }
      if (descriptions_.empty())  descriptions_.resize(entries.size());
      casadi_assert(descriptions_.size()==entries.size());
    }

    std::string IOSchemeCustomInternal::name() const {
      return "customIO";
    }

    std::string IOSchemeCustomInternal::entryNames() const {
      std::stringstream ss;
      for (int i=0;i<entries_.size();++i) {
         if (i!=0) ss << ", ";
         ss << entries_[i];
      }
      return ss.str();
    }

    std::string IOSchemeCustomInternal::entry(int i) const {
      casadi_assert_message(i>=0 && i<entries_.size(),
                            "customIO::entry(): requesting entry for index " << i
                            << ", but IOScheme is only length " << entries_.size());
      return entries_[i];
    }

    std::string IOSchemeCustomInternal::entryEnum(int i) const {
      return "";
    }

    std::string IOSchemeCustomInternal::describe(int i) const {
      if (descriptions_[i].empty()) {
        return entry(i);
      }  else {
        std::stringstream ss;
        ss << entry(i) <<  " '" << descriptions_[i] << "'";
        return ss.str();
      }
    }

    int IOSchemeCustomInternal::index(const std::string &name) const {
      std::map<std::string, int>::const_iterator it = entrymap_.find(name);
      casadi_assert_message(it!=entrymap_.end(), "customIO::index(): entry '" << name
                            << "' not available. Available entries are " << entryNames());
      return it->second;
    }

    int IOSchemeCustomInternal::size() const {
      return entries_.size();
    }

    void IOSchemeCustomInternal::print(std::ostream &stream) const {
      stream << "customIO(" << entryNames() << ")";
    }

    void IOSchemeCustomInternal::repr(std::ostream &stream) const {
      stream << "customIO(" << entryNames() << ")";
    }

} // namespace casadi
