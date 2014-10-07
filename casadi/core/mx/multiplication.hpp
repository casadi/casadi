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


#ifndef CASADI_MULTIPLICATION_HPP
#define CASADI_MULTIPLICATION_HPP

#include "mx_node.hpp"

/// \cond INTERNAL

namespace casadi {
  /** \brief An MX atomic for matrix-matrix product,
             note that the first factor must be provided transposed
      \author Joel Andersson
      \date 2010
  */
  template<bool TrX, bool TrY>
  class CASADI_CORE_EXPORT Multiplication : public MXNode {
  public:

    /** \brief  Constructor */
    Multiplication(const MX& z, const MX& x, const MX& y);

    /** \brief  Destructor */
    virtual ~Multiplication() {}

    /** \brief  Clone function */
    virtual Multiplication* clone() const { return new Multiplication(*this);}

    /** \brief  Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;

    /** \brief Generate code for the operation */
    virtual void generateOperation(std::ostream &stream, const std::vector<std::string>& arg,
                                   const std::vector<std::string>& res, CodeGenerator& gen) const;

    /// Evaluate the function (template)
    template<typename T, typename MatV, typename MatVV>
    void evaluateGen(const MatV& input, MatV& output, std::vector<int>& itmp, std::vector<T>& rtmp);

    /// Evaluate the function numerically
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp,
                           std::vector<double>& rtmp);

    /// Evaluate the function symbolically (SX)
    virtual void evaluateSX(const SXPtrV& input, SXPtrV& output, std::vector<int>& itmp,
                            std::vector<SXElement>& rtmp);

    /** \brief  Evaluate the function symbolically (MX) */
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed,
                            MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens,
                            bool output_given);

    /** \brief  Propagate sparsity */
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);

    /** \brief Get the operation */
    virtual int getOp() const { return OP_MATMUL;}

    /// Can the operation be performed inplace (i.e. overwrite the result)
    virtual int numInplace() const { return 1;}

    /** \brief Check if two nodes are equivalent up to a given depth */
    virtual bool isEqual(const MXNode* node, int depth) const
    { return sameOpAndDeps(node, depth) && dynamic_cast<const Multiplication<TrX, TrY>*>(node)!=0;}

    /// Helper class
    template<bool Tr>
    static MX tr(const MX& x) { return Tr ? x.T() : x;}
  };


  /** \brief An MX atomic for matrix-matrix product,
             note that the factor must be provided transposed
      \author Joel Andersson
      \date 2010
  */
  template<bool TrX, bool TrY>
  class CASADI_CORE_EXPORT DenseMultiplication : public Multiplication<TrX, TrY>{
  public:

    /** \brief  Constructor */
    DenseMultiplication(const MX& z, const MX& x, const MX& y)
        : Multiplication<TrX, TrY>(z, x, y) {}

    /** \brief  Destructor */
    virtual ~DenseMultiplication() {}

    /** \brief  Clone function */
    virtual DenseMultiplication* clone() const { return new DenseMultiplication(*this);}

    /** \brief Generate code for the operation */
    virtual void generateOperation(std::ostream &stream, const std::vector<std::string>& arg,
                                   const std::vector<std::string>& res, CodeGenerator& gen) const;
  };


} // namespace casadi
/// \endcond

#endif // CASADI_MULTIPLICATION_HPP
