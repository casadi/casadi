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

#ifndef MULTIPLICATION_HPP
#define MULTIPLICATION_HPP

#include "mx_node.hpp"

namespace CasADi{
/** \brief Represents matrix products
  \author Joel Andersson 
  \date 2010
  */
class Multiplication : public MXNode{
  public:
    
    /** \brief  Constructor */
    Multiplication(const MX& x, const MX& y);

    /** \brief  Destructor */
    virtual ~Multiplication(){}

    /** \brief  Clone function */
    virtual Multiplication* clone() const;

    /** \brief  Print */
    virtual void print(std::ostream &stream, const std::vector<std::string>& args) const;

    /** \brief  Evaluate the function and store the result in the node */
    virtual void evaluate(const std::vector<DMatrix*>& input, DMatrix& output, const VVDptr& fwdSeed, std::vector<DMatrix*>& fwdSens, const std::vector<DMatrix*>& adjSeed, VVDptr& adjSens, int nfwd, int nadj);

    /** \brief  Evaluate dense-dense */
    void evaluateDenseDense(const std::vector<DMatrix*>& input, DMatrix& output, const VVDptr& fwdSeed, std::vector<DMatrix*>& fwdSens, const std::vector<DMatrix*>& adjSeed, VVDptr& adjSens, int nfwd, int nadj);

    /** \brief  Evaluate sparse-sparse */
    void evaluateSparseSparse(const std::vector<DMatrix*>& input, DMatrix& output, const VVDptr& fwdSeed, std::vector<DMatrix*>& fwdSens, const std::vector<DMatrix*>& adjSeed, VVDptr& adjSens, int nfwd, int nadj);

    /** \brief  Evaluate symbolically (SX) */
    virtual void evaluateSX(const std::vector<SXMatrix*> &input, SXMatrix& output);

    // Check if the arguments are dense
    bool x_dense_, y_dense_;

    // Mapping corresponding to the transpose of y (no need to form the transpose explicitly)
    std::vector<int> y_trans_map_;
  
    // Create the sparsity pattern for the matrix-matrix product
    std::vector< std::vector< std::pair<int,int> > > prod_map_;
    


};

} // namespace CasADi


#endif // MULTIPLICATION_HPP
