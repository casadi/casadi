/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#ifndef NORM_HPP
#define NORM_HPP

#include "mx_node.hpp"

namespace CasADi{

/** \brief Matrix and vector norms
  This base class and the derived classes represent matrix and vector norms that are intended to be used when formulating convex 
  optimization problems. Note that they are not intended to be evaluated numerically or differentiated, instead the idea is that
  they should be eliminated from the computational graph during a reformulation (cf. CVX software).

  \author Joel Andersson 
  \date 2010-2011
*/
class Norm : public MXNode{
  public:

    /** \brief  Constructor */
    Norm(const MX& x);

    /** \brief  Destructor */
    virtual ~Norm(){}

    /** \brief  Evaluate the function numerically */
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens);

    /** \brief  Evaluate the function symbolically (SX) */
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens);

    /** \brief  Evaluate the function symbolically (MX) */
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given);

    /** \brief  Propagate sparsity */
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);

};

/** \brief Represents a 2-norm
  \author Joel Andersson 
  \date 2010
*/
class Norm2 : public Norm{
  public:

    /** \brief  Constructor */
    Norm2(const MX& x);

    /** \brief  Destructor */
    virtual ~Norm2(){}

    /** \brief  Clone function */
    virtual Norm2* clone() const;

    /** \brief  Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;

    /** \brief Get the operation */
    virtual int getOp() const{ return OP_NORM2;}
};

/** \brief Represents a Frobenius norm
  \author Joel Andersson 
  \date 2010
*/
class NormF : public Norm{
  public:

    /** \brief  Constructor */
    NormF(const MX& x);

    /** \brief  Destructor */
    virtual ~NormF(){}

    /** \brief  Clone function */
    virtual NormF* clone() const;

    /** \brief  Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;

    /** \brief Get the operation */
    virtual int getOp() const{ return OP_NORMF;}
};

/** \brief 1-norm
  \author Joel Andersson 
  \date 2010
*/
class Norm1 : public Norm{
  public:

    /** \brief  Constructor */
    Norm1(const MX& x);

    /** \brief  Destructor */
    virtual ~Norm1(){}

    /** \brief  Clone function */
    virtual Norm1* clone() const;

    /** \brief  Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;

    /** \brief Get the operation */
    virtual int getOp() const{ return OP_NORM1;}
};

/** \brief Represents an infinity-norm operation on a MX
  \author Joel Andersson 
  \date 2010
*/
class NormInf : public Norm{
  public:

    /** \brief  Constructor */
    NormInf(const MX& x);

    /** \brief  Destructor */
    virtual ~NormInf(){}

    /** \brief  Clone function */
    virtual NormInf* clone() const;

    /** \brief  Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;

    /** \brief Get the operation */
    virtual int getOp() const{ return OP_NORMINF;}
};

} // namespace CasADi

#endif // NORM_HPP
