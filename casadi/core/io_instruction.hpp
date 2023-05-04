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


#ifndef CASADI_IO_INSTRUCTION_HPP
#define CASADI_IO_INSTRUCTION_HPP

#include "mx_node.hpp"

namespace casadi {
  /** \brief An input or output instruction

      \author Joel Andersson
      \date 2017

      \identifier{1p6} */
  class CASADI_EXPORT IOInstruction : public MXNode {
  protected:
    // Input/output index
    casadi_int ind_;

    // Segment number
    casadi_int segment_;

    // Nonzero offset
    casadi_int offset_;

    // Constructor (called from derived classes)
    IOInstruction(casadi_int ind, casadi_int segment, casadi_int offset)
      : ind_(ind), segment_(segment), offset_(offset) {}

  public:
    /// Destructor
    ~IOInstruction() override {}

    // Get IO index
    casadi_int ind() const override { return ind_;}

    // Get IO segment
    casadi_int segment() const override { return segment_;}

    // Get IO offset
    casadi_int offset() const override { return offset_;}

    /** Obtain information about node */
    Dict info() const override;

    /** \brief Serialize an object without type information

        \identifier{1p7} */
    void serialize_body(SerializingStream& s) const override;

  protected:
    /** \brief Deserializing constructor

        \identifier{1p8} */
    explicit IOInstruction(DeserializingStream&s);
  };

  /** \brief Input instruction

      \identifier{1p9} */
  class CASADI_EXPORT Input : public IOInstruction {
  public:
    // Constructor (called from derived classes)
    Input(const Sparsity& sp, casadi_int ind, casadi_int segment, casadi_int offset);

    /// Destructor
    ~Input() override {}

    /** \brief Get the operation

        \identifier{1pa} */
    casadi_int op() const override { return OP_INPUT;}

    /** \brief  Print expression

        \identifier{1pb} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Generate code for the operation

        \identifier{1pc} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief Deserialize without type information

        \identifier{1pd} */
    static MXNode* deserialize(DeserializingStream& s) { return new Input(s); }

  protected:
    /** \brief Deserializing constructor

        \identifier{1pe} */
    explicit Input(DeserializingStream&s) : IOInstruction(s) {}
  };

  /** \brief Input instruction

      \identifier{1pf} */
  class CASADI_EXPORT Output : public IOInstruction {
  public:
    // Constructor (called from derived classes)
    Output(const MX& x, casadi_int ind, casadi_int segment, casadi_int offset);

    /// Destructor
    ~Output() override {}

    /** \brief  Number of outputs

        \identifier{1pg} */
    casadi_int nout() const override { return 0;}

    /** \brief Get the operation

        \identifier{1ph} */
    casadi_int op() const override { return OP_OUTPUT;}

    /** \brief  Print expression

        \identifier{1pi} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Generate code for the operation

        \identifier{1pj} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief Deserialize without type information

        \identifier{1pk} */
    static MXNode* deserialize(DeserializingStream& s) { return new Output(s); }
  protected:
    /** \brief Deserializing constructor

        \identifier{1pl} */
    explicit Output(DeserializingStream&s) : IOInstruction(s) {}
  };


  /// \endcond
} // namespace casadi

#endif // CASADI_IO_INSTRUCTION_HPP
