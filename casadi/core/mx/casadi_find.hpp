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


#ifndef CASADI_FIND_HPP
#define CASADI_FIND_HPP

#include "mx_node.hpp"
/// \cond INTERNAL

namespace casadi {
  /** \brief Finds the first nonzero element in a vector
      \author Joel Andersson
      \date 2015
  */
  class CASADI_EXPORT Find : public MXNode {
  public:
    /** \brief  Constructor */
    Find(const MX& x);

    /** \brief  Destructor */
    virtual ~Find() {}

    /** \brief  Print expression */
    virtual std::string print(const std::vector<std::string>& arg) const;

    /// Evaluate the function numerically
    virtual void eval(const double** arg, double** res, int* iw, double* w, int mem) const;

    /** \brief  Evaluate symbolically (MX) */
    virtual void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res);

    /** \brief Calculate forward mode directional derivatives */
    virtual void evalFwd(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void evalAdj(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens);

    /** \brief  Propagate sparsity forward */
    virtual void spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief Get the operation */
    virtual int op() const { return OP_FIND;}

    /** \brief Generate code for the operation */
    virtual void generate(CodeGenerator& g, const std::string& mem,
                          const std::vector<int>& arg, const std::vector<int>& res) const;
  };

} // namespace casadi

/// \endcond

#endif // CASADI_FIND_HPP
