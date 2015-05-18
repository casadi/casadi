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

    /** \brief  Print expression */
    virtual std::string print(const std::vector<std::string>& arg) const;

    /// Evaluate the function (template)
    template<typename T>
    void evalGen(const T* const* arg, T* const* res, int* itmp, T* rtmp);

    /// Evaluate the function numerically
    virtual void evalD(const double** input, double** output, int* itmp, double* rtmp);

    /// Evaluate the function symbolically (SX)
    virtual void evalSX(const SXElement** input, SXElement** output,
                            int* itmp, SXElement* rtmp);

    /** \brief  Evaluate symbolically (MX) */
    virtual void evalMX(const std::vector<MX>& arg, std::vector<MX>& res);

    /** \brief Calculate forward mode directional derivatives */
    virtual void evalFwd(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void evalAdj(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens);

    /** \brief Generate code for the operation */
    void generate(const std::vector<int>& arg, const std::vector<int>& res,
                  CodeGenerator& g) const;

    /** \brief  Propagate sparsity forward */
    virtual void spFwd(const bvec_t** arg, bvec_t** res, int* itmp, bvec_t* rtmp);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdj(bvec_t** arg, bvec_t** res, int* itmp, bvec_t* rtmp);

    /** \brief Get the operation */
    virtual int getOp() const { return OP_SET_SPARSE;}

    /** \brief Get number of temporary variables needed */
    virtual void nwork(size_t& n_arg, size_t& n_res, size_t& n_iw, size_t& n_w) const {
      n_arg=0;
      n_res=0;
      n_iw=0;
      n_w=size1();
    }
  };

} // namespace casadi

/// \endcond

#endif // CASADI_SET_SPARSE_HPP
