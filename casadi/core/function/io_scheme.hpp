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

  // Forward declaration
  class IOSchemeInternal;

  /** \brief Class with mapping between names and indices
   *
   * \author Joris Gillis
   * \date 2013
   */
class CASADI_CORE_EXPORT IOScheme : public SharedObject {

  public:

    /// Default constructor
    IOScheme();

    /// Constructor with enum
    IOScheme(InputOutputScheme scheme);

    /// Constructor with entry names
    IOScheme(const std::vector<std::string> &entries,
             const std::vector<std::string> &descriptions=std::vector<std::string>());

    #ifndef SWIGPYTHON
    IOScheme(
        const std::string &arg_s0, const std::string &arg_s1="",
        const std::string &arg_s2="", const std::string &arg_s3="",
        const std::string &arg_s4="", const std::string &arg_s5="",
        const std::string &arg_s6="", const std::string &arg_s7="",
        const std::string &arg_s8="", const std::string &arg_s9="",
        const std::string &arg_s10="", const std::string &arg_s11="",
        const std::string &arg_s12="", const std::string &arg_s13="",
        const std::string &arg_s14="", const std::string &arg_s15="",
        const std::string &arg_s16="", const std::string &arg_s17="",
        const std::string &arg_s18="", const std::string &arg_s19="");
   #endif // SWIG

    /** \brief  Access functions of the node */
    IOSchemeInternal* operator->();
    const IOSchemeInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);

    /// Name of the scheme
    std::string name() const;

    /// List available entries
    std::string entryNames() const;

    /// Get index by entry name
    int index(const std::string &name) const;

    /// Number of entries
    int size() const;

    /// Get the entry name by index
    std::string entry(int i) const;

    /** \brief Get the entry label by index
    * If scheme is unknown, returns the index as a string
    */
    std::string entryLabel(int i) const;

    /// Get the entry enum name by index
    std::string entryEnum(int i) const;

    /// Describe the index as an input
    std::string describeInput(int i) const;

    /// Describe the index as an output
    std::string describeOutput(int i) const;

    /// Describe the index
    std::string describe(int i) const;

    /// Check whether the scheme is known
    bool known() const;

    /// Check whether this scheme is compatible with the given size
    int compatibleSize(int size) const;

    #ifndef SWIG
    /// Print a description of the object
    void print(std::ostream &stream=std::cout) const;

    /// Print a representation of the object
    void repr(std::ostream &stream) const;
    #endif // SWIG

    #ifndef SWIGPYTHON
    template<class M>
    IOSchemeVector<M> operator()(
        const std::string arg_s0="", M arg_m0=M(),
        const std::string arg_s1="", M arg_m1=M(),
        const std::string arg_s2="", M arg_m2=M(),
        const std::string arg_s3="", M arg_m3=M(),
        const std::string arg_s4="", M arg_m4=M(),
        const std::string arg_s5="", M arg_m5=M(),
        const std::string arg_s6="", M arg_m6=M(),
        const std::string arg_s7="", M arg_m7=M(),
        const std::string arg_s8="", M arg_m8=M(),
        const std::string arg_s9="", M arg_m9=M(),
        const std::string arg_s10="", M arg_m10=M(),
        const std::string arg_s11="", M arg_m11=M(),
        const std::string arg_s12="", M arg_m12=M(),
        const std::string arg_s13="", M arg_m13=M(),
        const std::string arg_s14="", M arg_m14=M(),
        const std::string arg_s15="", M arg_m15=M(),
        const std::string arg_s16="", M arg_m16=M(),
        const std::string arg_s17="", M arg_m17=M(),
        const std::string arg_s18="", M arg_m18=M(),
        const std::string arg_s19="", M arg_m19=M()) {
      std::vector<std::string> k;
      std::vector<M> v;
      if (arg_s0!="") { k.push_back(arg_s0);  v.push_back(arg_m0); }
      if (arg_s1!="") { k.push_back(arg_s1);  v.push_back(arg_m1); }
      if (arg_s2!="") { k.push_back(arg_s2);  v.push_back(arg_m2); }
      if (arg_s3!="") { k.push_back(arg_s3);  v.push_back(arg_m3); }
      if (arg_s4!="") { k.push_back(arg_s4);  v.push_back(arg_m4); }
      if (arg_s5!="") { k.push_back(arg_s5);  v.push_back(arg_m5); }
      if (arg_s6!="") { k.push_back(arg_s6);  v.push_back(arg_m6); }
      if (arg_s7!="") { k.push_back(arg_s7);  v.push_back(arg_m7); }
      if (arg_s8!="") { k.push_back(arg_s8);  v.push_back(arg_m8); }
      if (arg_s9!="") { k.push_back(arg_s9);  v.push_back(arg_m9); }
      if (arg_s10!="") { k.push_back(arg_s10);  v.push_back(arg_m10); }
      if (arg_s11!="") { k.push_back(arg_s11);  v.push_back(arg_m11); }
      if (arg_s12!="") { k.push_back(arg_s12);  v.push_back(arg_m12); }
      if (arg_s13!="") { k.push_back(arg_s13);  v.push_back(arg_m13); }
      if (arg_s14!="") { k.push_back(arg_s14);  v.push_back(arg_m14); }
      if (arg_s15!="") { k.push_back(arg_s15);  v.push_back(arg_m15); }
      if (arg_s16!="") { k.push_back(arg_s16);  v.push_back(arg_m16); }
      if (arg_s17!="") { k.push_back(arg_s17);  v.push_back(arg_m17); }
      if (arg_s18!="") { k.push_back(arg_s18);  v.push_back(arg_m18); }
      if (arg_s19!="") { k.push_back(arg_s19);  v.push_back(arg_m19); }
      return operator()(k, v);
    }
    #endif // SWIG

  #ifndef SWIGPYTHON
  template<class M>
  std::vector<M> operator()(
      const std::vector<M> arg_m0,
      const std::string &arg_s0="", const std::string &arg_s1="",
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
    if (arg_s4!="") { k.push_back(arg_s4); }
    if (arg_s5!="") { k.push_back(arg_s5); }
    if (arg_s6!="") { k.push_back(arg_s6); }
    if (arg_s7!="") { k.push_back(arg_s7); }
    if (arg_s8!="") { k.push_back(arg_s8); }
    if (arg_s9!="") { k.push_back(arg_s9); }
    if (arg_s10!="") { k.push_back(arg_s10); }
    if (arg_s11!="") { k.push_back(arg_s11); }
    if (arg_s12!="") { k.push_back(arg_s12); }
    if (arg_s13!="") { k.push_back(arg_s13); }
    if (arg_s14!="") { k.push_back(arg_s14);  }
    if (arg_s15!="") { k.push_back(arg_s15);  }
    if (arg_s16!="") { k.push_back(arg_s16);  }
    if (arg_s17!="") { k.push_back(arg_s17);  }
    if (arg_s18!="") { k.push_back(arg_s18);  }
    if (arg_s19!="") { k.push_back(arg_s19);   }
    return operator()(arg_m0, k);
  }
  #endif // SWIG

  template<class M>
    std::vector<M> operator()(const std::vector<M> arg_m,
                              const std::vector<std::string> & arg_s) {
    casadi_assert(arg_m.size()==size());
    std::vector<M> ret;
    for (int i=0;i<arg_s.size();++i) {
      int k = index(arg_s[i]);
      ret.push_back(arg_m[k]);
    }

    return ret;
  }

  template<class M>
    IOSchemeVector<M> operator()(const std::vector<std::string> arg_s,
                                 const std::vector<M> & arg_m) {
    casadi_assert(arg_s.size()==arg_m.size());
    std::vector<M> v(size());
    for (int i=0;i<arg_s.size();++i) {
      int k = index(arg_s[i]);
      v[k]=arg_m[i];
    }

    return IOSchemeVector<M>(v, *this);
  }

  template<class M>
    IOSchemeVector<M> fromVector(const std::vector<M> & arg_m) {
    casadi_assert(!known() || size()==arg_m.size());
    return IOSchemeVector<M>(arg_m, *this);
  }
};

} // namespace casadi

#endif // CASADI_IO_SCHEME_HPP
