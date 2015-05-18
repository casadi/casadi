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


#ifndef CASADI_TRANSPOSE_HPP
#define CASADI_TRANSPOSE_HPP

#include "mx_node.hpp"
#include <map>
#include <stack>

/// \cond INTERNAL

namespace casadi {
  /** \brief Matrix transpose
      \author Joel Andersson
      \date 2013
  */
  class CASADI_EXPORT Transpose : public MXNode {
  public:

    /// Constructor
    Transpose(const MX& x);

    /// Clone function
    virtual Transpose* clone() const { return new Transpose(*this);}

    /// Destructor
    virtual ~Transpose() {}

    /// Evaluate the function (template)
    template<typename T>
    void evalGen(const T* const* arg, T* const* res, int* iw, T* w);

    /// Evaluate the function numerically
    virtual void evalD(const double** arg, double** res, int* iw, double* w);

    /// Evaluate the function symbolically (SX)
    virtual void evalSX(const SXElement** arg, SXElement** res, int* iw, SXElement* w);

    /** \brief  Evaluate symbolically (MX) */
    virtual void evalMX(const std::vector<MX>& arg, std::vector<MX>& res);

    /** \brief Calculate forward mode directional derivatives */
    virtual void evalFwd(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void evalAdj(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens);

    /** \brief  Propagate sparsity forward */
    virtual void spFwd(const bvec_t** arg,
                       bvec_t** res, int* iw, bvec_t* w);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdj(bvec_t** arg,
                       bvec_t** res, int* iw, bvec_t* w);

    /** \brief  Print expression */
    virtual std::string print(const std::vector<std::string>& arg) const;

    /** \brief Generate code for the operation */
    virtual void generate(const std::vector<int>& arg, const std::vector<int>& res,
                          CodeGenerator& g) const;

    /** \brief Get the operation */
    virtual int getOp() const { return OP_TRANSPOSE;}

    /** \brief Get required length of iw field */
    virtual size_t sz_iw() const { return size2()+1;}

    /// Transpose
    virtual MX getTranspose() const { return dep();}

    /// Solve for square linear system
    //virtual MX getSolve(const MX& r, bool tr, const LinearSolver& linear_solver) const {
    // return dep()->getSolve(r, !tr, linear_solver);} // FIXME #1001

    /** \brief Check if two nodes are equivalent up to a given depth */
    virtual bool zz_isEqual(const MXNode* node, int depth) const {
      return sameOpAndDeps(node, depth);
    }
  };

  /** \brief Matrix transpose (dense)
      \author Joel Andersson
      \date 2013
  */
  class CASADI_EXPORT DenseTranspose : public Transpose {
  public:

    /// Constructor
    DenseTranspose(const MX& x) : Transpose(x) {}

    /// Clone function
    virtual DenseTranspose* clone() const { return new DenseTranspose(*this);}

    /// Destructor
    virtual ~DenseTranspose() {}

    /// Evaluate the function (template)
    template<typename T>
    void evalGen(const T* const* arg, T* const* res, int* iw, T* w);

    /// Evaluate the function numerically
    virtual void evalD(const double** arg, double** res, int* iw, double* w);

    /// Evaluate the function symbolically (SX)
    virtual void evalSX(const SXElement** arg, SXElement** res,
                        int* iw, SXElement* w);

    /** \brief  Propagate sparsity forward */
    virtual void spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w);

    /** \brief Generate code for the operation */
    virtual void generate(const std::vector<int>& arg, const std::vector<int>& res,
                          CodeGenerator& g) const;

    /** \brief Get required length of iw field */
    virtual size_t sz_iw() const { return 0;}
  };



} // namespace casadi

/// \endcond

#endif // CASADI_TRANSPOSE_HPP
