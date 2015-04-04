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


#ifndef CASADI_SET_SPARSE_HPP
#define CASADI_SET_SPARSE_HPP

#include "mx_node.hpp"
/// \cond INTERNAL

namespace casadi {
  /** \brief Change the sparsity of an expression
      \author Joel Andersson
      \date 2011-2013
  */
  class CASADI_EXPORT SetSparse : public MXNode {
  public:

    /** \brief  Constructor */
    SetSparse(const MX& x, const Sparsity& sp);

    /** \brief  Destructor */
    virtual ~SetSparse() {}

    /** \brief  Clone function */
    virtual SetSparse * clone() const;

    /** \brief  Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;

    /// Evaluate the function (template)
    template<typename T>
    void evalGen(const T* const* arg, T* const* res, int* itmp, T* rtmp);

    /// Evaluate the function numerically
    virtual void evalD(cp_double* input, p_double* output, int* itmp, double* rtmp);

    /// Evaluate the function symbolically (SX)
    virtual void evalSX(cp_SXElement* input, p_SXElement* output,
                            int* itmp, SXElement* rtmp);

    /** \brief  Evaluate the function symbolically (MX) */
    virtual void eval(const cpv_MX& input, const pv_MX& output);

    /** \brief Calculate forward mode directional derivatives */
    virtual void evalFwd(const std::vector<cpv_MX>& fwdSeed, const std::vector<pv_MX>& fwdSens);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void evalAdj(const std::vector<pv_MX>& adjSeed, const std::vector<pv_MX>& adjSens);

    /** \brief Generate code for the operation */
    void generate(std::ostream &stream, const std::vector<int>& arg,
                           const std::vector<int>& res, CodeGenerator& gen) const;

    /** \brief  Propagate sparsity forward */
    virtual void spFwd(cp_bvec_t* arg,
                       p_bvec_t* res, int* itmp, bvec_t* rtmp);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdj(p_bvec_t* arg,
                       p_bvec_t* res, int* itmp, bvec_t* rtmp);

    /** \brief Get the operation */
    virtual int getOp() const { return OP_SET_SPARSE;}

    /// Get number of temporary variables needed
    virtual void nTmp(size_t& ni, size_t& nr) { ni=0; nr=size1();}
  };

} // namespace casadi

/// \endcond

#endif // CASADI_SET_SPARSE_HPP
