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


#ifndef CASADI_REPMAT_HPP
#define CASADI_REPMAT_HPP

#include "mx_node.hpp"
#include <map>
#include <stack>

/// \cond INTERNAL

namespace casadi {

  /** \brief Horizontal repmat
      \author Joris Gillis
      \date 2015
  */
  class CASADI_EXPORT HorzRepmat : public MXNode {
  public:

    /// Constructor
    HorzRepmat(const MX& x, int n);

    /// Clone function
    virtual HorzRepmat* clone() const { return new HorzRepmat(*this);}

    /// Evaluate the function (template)
    template<typename T>
    void evalGen(const T** arg, T** res, int* iw, T* w);

    /// Destructor
    virtual ~HorzRepmat() {}

    /** \brief  Print expression */
    virtual std::string print(const std::vector<std::string>& arg) const;

    /// Evaluate the function numerically
    virtual void evalD(const double** arg, double** res, int* iw, double* w);

    /// Evaluate the function symbolically (SX)
    virtual void evalSX(const SXElement** arg, SXElement** res, int* iw, SXElement* w);

    /** \brief  Evaluate symbolically (MX) */
    virtual void evalMX(const std::vector<MX>& arg, std::vector<MX>& res);

    /** \brief  Propagate sparsity forward */
    virtual void spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w);

    /** \brief Calculate forward mode directional derivatives */
    virtual void evalFwd(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void evalAdj(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens);

    /** \brief Generate code for the operation */
    void generate(const std::vector<int>& arg, const std::vector<int>& res,
                  CodeGenerator& g) const;

    /** \brief Get the operation */
    virtual int getOp() const { return OP_HORZREPMAT;}

    int n_;
  };

  /** \brief Horizontal repsum
      \author Joris Gillis
      \date 2015
  */
  class CASADI_EXPORT HorzRepsum : public MXNode {
  public:

    /// Constructor
    HorzRepsum(const MX& x, int n);

    /// Clone function
    virtual HorzRepsum* clone() const { return new HorzRepsum(*this);}

    /// Evaluate the function (template)
    template<typename T, typename R>
    void evalGen(const T** arg, T** res, int* iw, T* w,
      R reduction);

    /// Destructor
    virtual ~HorzRepsum() {}

    /** \brief  Print expression */
    virtual std::string print(const std::vector<std::string>& arg) const;

    /// Evaluate the function numerically
    virtual void evalD(const double** arg, double** res, int* iw, double* w);

    /// Evaluate the function symbolically (SX)
    virtual void evalSX(const SXElement** arg, SXElement** res, int* iw, SXElement* w);

    /** \brief  Evaluate symbolically (MX) */
    virtual void evalMX(const std::vector<MX>& arg, std::vector<MX>& res);

    /** \brief  Propagate sparsity forward */
    virtual void spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w);

    /** \brief Calculate forward mode directional derivatives */
    virtual void evalFwd(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void evalAdj(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens);

    /** \brief Generate code for the operation */
    void generate(const std::vector<int>& arg, const std::vector<int>& res,
                  CodeGenerator& g) const;

    /** \brief Get the operation */
    virtual int getOp() const { return OP_HORZREPSUM;}

    int n_;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_REPMAT_HPP
