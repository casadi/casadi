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


#ifndef CASADI_OUTPUT_SX_HPP
#define CASADI_OUTPUT_SX_HPP

#include "sx_node.hpp"
#include "generic_shared_internal.hpp"
#include "generic_shared.hpp"

/// \cond INTERNAL
namespace casadi {

  class OutputSX;

  class CASADI_EXPORT SharedSXElem :
      public GenericShared<SharedSXElem, OutputSX> {
    public:

    static bool test_cast(const OutputSX* ptr) {
      return ptr;
    }

    using internal_base_type = OutputSX;
    using base_type = SharedSXElem;

  };

  class CASADI_EXPORT WeakRefSXElem :
      public GenericWeakRef<SharedSXElem, OutputSX> {
  public:
    WeakRefSXElem(int dummy=0) : GenericWeakRef<SharedSXElem, OutputSX>(dummy) {
    }
    WeakRefSXElem(SharedSXElem shared) : GenericWeakRef<SharedSXElem, OutputSX>(shared) {
    }
  };

  typedef GenericWeakRefInternal<SharedSXElem, OutputSX> WeakRefInternalSXElem;

    class CASADI_EXPORT OutputSX :
      public SXNode,
      public GenericSharedInternal<SharedSXElem, OutputSX> {
    friend class GenericShared<SharedSXElem, OutputSX>;
    friend class SharedObject;
    friend class GenericWeakRef<SharedSXElem, OutputSX>;
    friend class GenericSharedInternal<SharedSXElem, OutputSX>;
    friend class Memory;
    friend class UniversalNodeOwner;
  public:

    using weak_ref_type = WeakRefInternalSXElem;

    /** \brief  Constructor

        \identifier{294} */
    OutputSX(const SXElem& dep, int oind) : dep_(dep), oind_(oind) {
    }

    /// Empty constructor
    OutputSX() {

    }

    // Class name
    std::string class_name() const override {return "OutputSX";}

    /** \brief  Print expression

        \identifier{296} */
    std::string print(const std::string& arg1, const std::string& arg2) const override {
      return arg1 + "{" + str(oind_) + "}";
    }

    /** \brief  Destructor

        \identifier{297} */
    ~OutputSX() override {
      safe_delete(dep_.assignNoDelete(casadi_limits<SXElem>::nan));
    }

    /** \brief Get the operation

        \identifier{298} */
    casadi_int op() const override { return -1;}

    /** \brief  Number of dependencies

        \identifier{299} */
    casadi_int n_dep() const override { return 1;}

    /** \brief  get the reference of a dependency

        \identifier{29a} */
    const SXElem& dep(casadi_int i) const override { return dep_; }
    SXElem& dep(casadi_int i) override { return dep_; }

    /** \brief  The dependencies of the node

        \identifier{29b} */
    SXElem dep_;

    /** \brief  Output index

        \identifier{29c} */
    int oind_;

    /** \brief  Is the node an output node?

        \identifier{29d} */
    bool is_output() const override { return true; }

    /** \brief Get the index of evaluation output

        \identifier{29e} */
    casadi_int which_output() const override { return oind_; }

    static std::vector<SXElem> split(const SXElem& e, casadi_int n) {
      std::vector<SXElem> ret(n);
      for (casadi_int i=0;i<n;++i) {
        ret[i] = e.get_output(i);
      }
      return ret;
    }

    void serialize_node(SerializingStream& s) const override {
      s.pack("OutputSX::dep", dep_);
      s.pack("OutputSX::oind", oind_);
    }

    static SXNode* deserialize(DeserializingStream& s) {
      SXElem dep;
      int oind;
      s.unpack("OutputSX::dep", dep);
      s.unpack("OutputSX::oind", oind);
      return new OutputSX(dep, oind);
    }

  };


} // namespace casadi
/// \endcond

#endif // CASADI_OUTPUT_SX_HPP
