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

/** \brief Represents any type of general norm
  \author Joel Andersson 
  \date 2010
*/
class Norm : public MXNode{
public:

/** \brief  Constructor */
Norm(const MX& x);

/** \brief  Evaluate */
virtual void evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens, int nfwd);

/// Symbolic forward sensitivities
virtual MX adFwd(const std::vector<MX>& jx);

};

/** \brief Represents a 2-norm operation on a MX
  Frobenius norm
  \author Joel Andersson 
  \date 2010
*/
class Norm2 : public Norm{
public:

/** \brief  Constructor */
Norm2(const MX& x);

/** \brief  Clone function */
virtual Norm2* clone() const;

/** \brief  Print */
virtual void print(std::ostream &stream, const std::vector<std::string>& args) const;

/** \brief  Evaluate */
virtual void evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens, int nfwd);

/** \brief Symbolic forward sensitivities.  */
virtual MX adFwd(const std::vector< MX > & jx	);

};

/** \brief Represents a 2-norm squared operation on a MX
  Frobenius norm
  \author Joel Andersson 
  \date 2010
*/
class Norm22 : public Norm{
public:

/** \brief  Constructor */
Norm22(const MX& x);

/** \brief  Clone function */
virtual Norm22* clone() const;

/** \brief  Print */
virtual void print(std::ostream &stream, const std::vector<std::string>& args) const;

/** \brief  Evaluate */
virtual void evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens, int nfwd);

/** \brief Symbolic forward sensitivities. */
virtual MX adFwd(const std::vector< MX > & jx	);

};

/** \brief Represents a 1-norm operation on a MX
  Entrywise Norm
  \author Joel Andersson 
  \date 2010
*/
class Norm1 : public Norm{
public:

/** \brief  Constructor */
Norm1(const MX& x);

/** \brief  Clone function */
virtual Norm1* clone() const;

/** \brief  Print */
virtual void print(std::ostream &stream, const std::vector<std::string>& args) const;

/** \brief  Evaluate */
virtual void evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens, int nfwd);

};

/** \brief Represents an infinity-norm operation on a MX
  \author Joel Andersson 
  \date 2010
*/
class NormInf : public Norm{
public:

/** \brief  Constructor */
NormInf(const MX& x);

/** \brief  Clone function */
virtual NormInf* clone() const;

/** \brief  Print */
virtual void print(std::ostream &stream, const std::vector<std::string>& args) const;

/** \brief  Evaluate */
virtual void evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens, int nfwd);


};

} // namespace CasADi

#endif // NORM_HPP
