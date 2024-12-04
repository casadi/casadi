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

/// \cond INTERNAL
namespace casadi {

    class CASADI_EXPORT OutputSX : public SXNode {
  private:

    /** \brief  Constructor */
    OutputSX(const SXElem& dep, int oind) : dep_(dep), oind_(oind) {

    }

  public:
    /** \brief  Create a unary expression */
    inline static SXElem create(const SXElem& dep, int oind) {
      return SXElem::create(new OutputSX(dep, oind));
    }

    // Class name
    std::string class_name() const override {return "OutputSX";}

    /** \brief  Print expression */
    std::string print(const std::string& arg1, const std::string& arg2) const override {
      return "output(" + arg1 + "," + str(oind_) + ")";
    }

    /** \brief  Destructor */
    ~OutputSX() override {
      safe_delete(dep_.assignNoDelete(casadi_limits<SXElem>::nan));
    }

    /** \brief Get the operation */
    casadi_int op() const override { return -1;}

    /** \brief  Number of dependencies */
    casadi_int n_dep() const override { return 1;}

    /** \brief  get the reference of a dependency */
    const SXElem& dep(casadi_int i) const override { return dep_; }
    SXElem& dep(casadi_int i) override { return dep_; }

    /** \brief  The dependencies of the node */
    SXElem dep_;

    /** \brief  Output index */
    int oind_;

    static std::vector<SXElem> split(const SXElem& e, casadi_int n) {
      std::vector<SXElem> ret(n);
      for (casadi_int i=0;i<n;++i) {
        ret[i] = OutputSX::create(e, i);
        // Optimization: get_output() cached like in MultipleOutput
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
