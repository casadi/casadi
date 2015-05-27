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


#ifndef CASADI_IO_SCHEME_HPP
#define CASADI_IO_SCHEME_HPP

#include "../shared_object.hpp"
#include "schemes_metadata.hpp"
#include <string>
#include <vector>
#include <iostream>

namespace casadi {
  inline std::vector<std::string>
  IOScheme(InputOutputScheme scheme) {
    int sz = getSchemeSize(scheme);
    std::vector<std::string> data(sz);
    for (int i=0; i<sz; ++i) {
      data[i] = getSchemeEntryName(scheme, i) + ":"
        + getSchemeEntryEnumName(scheme, i)
        + " aka '" + getSchemeEntryName(scheme, i) + "'";
    }
    return data;
  }

  inline std::vector<std::string>
  IOScheme(const std::vector<std::string> &entries=std::vector<std::string>(),
           const std::vector<std::string> &descriptions=std::vector<std::string>()) {
    if (descriptions.empty()) {
      return entries;
    } else {
      casadi_assert(entries.size()==descriptions.size());
      std::vector<std::string> data=entries;
      for (size_t i=0; i!=entries.size(); ++i) {
        data[i] += ":" + descriptions[i];
      }
      return data;
    }
  }

#ifndef SWIGPYTHON
  inline std::vector<std::string>
  IOScheme(const std::string &arg_s0, const std::string &arg_s1="",
           const std::string &arg_s2="", const std::string &arg_s3="",
           const std::string &arg_s4="", const std::string &arg_s5="",
           const std::string &arg_s6="", const std::string &arg_s7="",
           const std::string &arg_s8="", const std::string &arg_s9="",
           const std::string &arg_s10="", const std::string &arg_s11="",
           const std::string &arg_s12="", const std::string &arg_s13="",
           const std::string &arg_s14="", const std::string &arg_s15="",
           const std::string &arg_s16="", const std::string &arg_s17="",
           const std::string &arg_s18="", const std::string &arg_s19="") {
    std::vector<std::string> k;
    if (arg_s0!="") { k.push_back(arg_s0);}
    if (arg_s1!="") { k.push_back(arg_s1);}
    if (arg_s2!="") { k.push_back(arg_s2);}
    if (arg_s3!="") { k.push_back(arg_s3);}
    if (arg_s4!="") { k.push_back(arg_s4);}
    if (arg_s5!="") { k.push_back(arg_s5);}
    if (arg_s6!="") { k.push_back(arg_s6);}
    if (arg_s7!="") { k.push_back(arg_s7);}
    if (arg_s8!="") { k.push_back(arg_s8);}
    if (arg_s9!="") { k.push_back(arg_s9);}
    if (arg_s10!="") { k.push_back(arg_s10);}
    if (arg_s11!="") { k.push_back(arg_s11);}
    if (arg_s12!="") { k.push_back(arg_s12);}
    if (arg_s13!="") { k.push_back(arg_s13);}
    if (arg_s14!="") { k.push_back(arg_s14);}
    if (arg_s15!="") { k.push_back(arg_s15);}
    if (arg_s16!="") { k.push_back(arg_s16);}
    if (arg_s17!="") { k.push_back(arg_s17);}
    if (arg_s18!="") { k.push_back(arg_s18);}
    if (arg_s19!="") { k.push_back(arg_s19);}
    return k;
  }
#endif // SWIG

} // namespace casadi

#endif // CASADI_IO_SCHEME_HPP
