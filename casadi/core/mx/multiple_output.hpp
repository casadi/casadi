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


#ifndef CASADI_MULTIPLE_OUTPUT_HPP
#define CASADI_MULTIPLE_OUTPUT_HPP

#include "mx_node.hpp"
#include "../function/function.hpp"
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

    /** \brief  Constructor */
    MultipleOutput();

    /** \brief  Destructor */
    virtual ~MultipleOutput();

    /** \brief  Number of outputs */
    virtual int nout() const=0;

    /** \brief  Get an output */
    virtual MX getOutput(int oind) const;

    /** \brief  Get the sparsity of output oind */
    virtual const Sparsity& sparsity(int oind) const=0;

    /** \brief  Check if a multiple output node */
    virtual bool isMultipleOutput() const {return true;}

  };

  class CASADI_EXPORT OutputNode : public MXNode {
  public:

    /** \brief  Constructor */
    OutputNode(const MX& parent, int oind);

    /** \brief  Destructor */
    virtual ~OutputNode();

    /** \brief  Print expression */
    virtual std::string print(const std::vector<std::string>& arg) const;

    /** \brief Is the node nonlinear */
    virtual bool isNonLinear() {return true;}

    /** \brief  Check if evaluation output */
    virtual bool isOutputNode() const {return true;}

    /** \brief  Get function input */
    virtual int getFunction_input() const { return -1;}

    /** \brief  Get function output */
    virtual int getFunctionOutput() const { return oind_;}

    /** \brief Get the operation */
    virtual int op() const { return -1;}

    /// Create a horizontal concatenation node
    virtual MX getHorzcat(const std::vector<MX>& x) const { return dep()->getHorzcat(x);}

    /// Create a vertical concatenation node (vectors only)
    virtual MX getVertcat(const std::vector<MX>& x) const { return dep()->getVertcat(x);}

    /** \brief  Output index */
    int oind_;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_MULTIPLE_OUTPUT_HPP
