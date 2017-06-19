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


#ifndef CASADI_IO_INSTRUCTION_HPP
#define CASADI_IO_INSTRUCTION_HPP

#include "mx_node.hpp"

namespace casadi {
  /** \brief An input or output instruction
      \author Joel Andersson
      \date 2017
  */
  class CASADI_EXPORT IOInstruction : public MXNode {
  protected:
    // Input/output index
    int ind_;

    // Segment number
    int segment_;

    // Nonzero offset
    int offset_;

    // Constructor (called from derived classes)
    IOInstruction(int ind, int segment, int offset)
      : ind_(ind), segment_(segment), offset_(offset) {}

  public:
    /// Destructor
    ~IOInstruction() override {}

    // Get IO index
    int ind() const override { return ind_;}

    // Get IO segment
    int segment() const override { return segment_;}

    // Get IO offset
    int offset() const override { return offset_;}
  };

  /** \brief Input instruction  */
  class CASADI_EXPORT Input : public IOInstruction {
  public:
    // Constructor (called from derived classes)
    Input(const Sparsity& sp, int ind, int segment, int offset);

    /// Destructor
    ~Input() override {}

    /** \brief Get the operation */
    int op() const override { return OP_INPUT;}

    /** \brief  Print expression */
    std::string print(const std::vector<std::string>& arg) const override;

    /** \brief Generate code for the operation */
    void generate(CodeGenerator& g, const std::string& mem,
                  const std::vector<int>& arg, const std::vector<int>& res) const override;
  };

  /** \brief Input instruction  */
  class CASADI_EXPORT Output : public IOInstruction {
  public:
    // Constructor (called from derived classes)
    Output(const MX& x, int ind, int segment, int offset);

    /// Destructor
    ~Output() override {}

    /** \brief  Number of outputs */
    int nout() const override { return 0;}

    /** \brief Get the operation */
    int op() const override { return OP_OUTPUT;}

    /** \brief  Print expression */
    std::string print(const std::vector<std::string>& arg) const override;

    /** \brief Generate code for the operation */
    void generate(CodeGenerator& g, const std::string& mem,
                  const std::vector<int>& arg, const std::vector<int>& res) const override;
  };


  /// \endcond
} // namespace casadi

#endif // CASADI_IO_INSTRUCTION_HPP
