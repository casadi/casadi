/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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


#ifndef CASADI_MULTIPLE_OUTPUT_HPP
#define CASADI_MULTIPLE_OUTPUT_HPP

#include "mx_node.hpp"
#include "function.hpp"
#include <set>

/// \cond INTERNAL

namespace casadi {

  /// Forward declaration
  class OutputNode;

  /**
      \author Joel Andersson
      \date 2010
  */
  class CASADI_EXPORT MultipleOutput : public MXNode {
    friend class OutputNode;
  public:

    /** \brief  Constructor

        \identifier{1pm} */
    MultipleOutput();

    /** \brief  Destructor

        \identifier{1pn} */
    ~MultipleOutput() override;

    /** \brief  Number of outputs

        \identifier{1po} */
    casadi_int nout() const override=0;

    /** \brief  Get an output

        \identifier{1pp} */
    MX get_output(casadi_int oind) const override;

    /** \brief  Get the sparsity of output oind

        \identifier{1pq} */
    const Sparsity& sparsity(casadi_int oind) const override=0;

    /** \brief  Check if a multiple output node

        \identifier{1pr} */
    bool has_output() const override {return true;}

  protected:
    /** \brief Deserializing constructor

        \identifier{1ps} */
    explicit MultipleOutput(DeserializingStream& s) : MXNode(s) {}

  };

  class CASADI_EXPORT OutputNode : public MXNode {
  public:

    /** \brief  Constructor

        \identifier{1pt} */
    OutputNode(const MX& parent, casadi_int oind);

    /** \brief  Destructor

        \identifier{1pu} */
    ~OutputNode() override;

    /** \brief  Print expression

        \identifier{1pv} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief  Check if evaluation output

        \identifier{1pw} */
    bool is_output() const override {return true;}

    /** \brief  Get function output

        \identifier{1px} */
    casadi_int which_output() const override { return oind_;}

    /** \brief Get the operation

        \identifier{1py} */
    casadi_int op() const override { return -1;}

    /// Create a horizontal concatenation node
    MX get_horzcat(const std::vector<MX>& x) const override { return dep()->get_horzcat(x);}

    /// Create a vertical concatenation node (vectors only)
    MX get_vertcat(const std::vector<MX>& x) const override { return dep()->get_vertcat(x);}

    /** Obtain information about node */
    Dict info() const override { return {{"oind", oind_}}; }

    /** \brief  Output index

        \identifier{1pz} */
    casadi_int oind_;

    /** \brief Serialize an object without type information

        \identifier{1q0} */
    void serialize_body(SerializingStream& s) const override;

    /** \brief Deserialize without type information

        \identifier{1q1} */
    static MXNode* deserialize(DeserializingStream& s) { return new OutputNode(s); }

  protected:
    /** \brief Deserializing constructor

        \identifier{1q2} */
    explicit OutputNode(DeserializingStream& s);
  };

} // namespace casadi
/// \endcond

#endif // CASADI_MULTIPLE_OUTPUT_HPP
