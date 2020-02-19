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


#ifndef CASADI_BSPLINE_HPP
#define CASADI_BSPLINE_HPP

#include "mx_node.hpp"
#include <map>
#include <stack>

/// \cond INTERNAL

namespace casadi {
  /** \brief BSpline Node
      \author Joris Gillis
      \date 2017-2019
  */
  class CASADI_EXPORT BSpline : public MXNode {
  public:

    static MX create(const MX& x, const std::vector< MX >& knots,
          const MX& coeffs,
          const std::vector<casadi_int>& degree,
          casadi_int m,
          const Dict& opts);

    /// Constructor
    BSpline(const MX& x, const MX& knots,
            const std::vector<casadi_int>& offset,
            const MX& coeffs,
            const std::vector<casadi_int>& degree,
            casadi_int m,
            const std::vector<casadi_int>& lookup_mode);

    /// Destructor
    ~BSpline() override {}

    /** \brief  Print expression */
    std::string disp(const std::vector<std::string>& arg) const override;

    static void prepare(casadi_int m, const std::vector<casadi_int>& offset,
      const std::vector<casadi_int>& degree, casadi_int &coeffs_size,
      std::vector<casadi_int>& coeffs_dims, std::vector<casadi_int>& strides);

    static casadi_int get_coeff_size(casadi_int m, const std::vector<casadi_int>& offset,
      const std::vector<casadi_int>& degree);

    MX knots_;
    // Numeric coefficients
    MX coeffs_;
    std::vector<double> knots_vec_;
    std::vector<double> coeffs_vec_;
    const std::vector<double>* knots_ptr_;
    const std::vector<double>* coeffs_ptr_;

    std::vector<casadi_int> offset_;
    std::vector<casadi_int> degree_;
    casadi_int m_;
    std::vector<casadi_int> lookup_mode_;

    // Derived fiels
    std::vector<casadi_int> strides_;
    std::vector<casadi_int> coeffs_dims_;
    casadi_int coeffs_size_;

 

    /** \brief Jacobian
     * 
     * Derivatives are computed by transforming the coefficient matrix
     * This is efficient
     */
    mutable MX jac_cache_;

    MX jac_cached() const;

    /** \brief Get required length of iw field */
    static size_t n_iw(const std::vector<casadi_int> &degree);

    /** \brief Get required length of w field */
    static size_t n_w(const std::vector<casadi_int> &degree);

    /** \brief Get required length of iw field */
    size_t sz_iw() const override;

    /** \brief Get required length of w field */
    size_t sz_w() const override;

    /** \brief Get the operation */
    casadi_int op() const override { return OP_BSPLINE;}

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /** \brief  Evaluate symbolically (MX) */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief Generate code for the operation */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief Generate code for the operation */
    std::string generate_coeff(CodeGenerator& g) const;

    /** \brief Generate code for the operation */
    std::string generate_grid(CodeGenerator& g) const;

    /** \brief Deserialize without type information */
    static MXNode* deserialize(DeserializingStream& s) { return new BSpline(s); }

    MX derivative_coeff(casadi_int i, const MX& coeffs, const MX& knots) const;

    MX jac(const MX& x, const MX& coeffs, const MX& kn) const;

    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream& s) const override;

    static MX dual(const MX& x,
          const std::vector< MX >& knots,
          const std::vector<casadi_int>& degree,
          const Dict& opts);


  protected:

    /** \brief Deserializing constructor */
    explicit BSpline(DeserializingStream& s);

  };


} // namespace casadi
/// \endcond

#endif // CASADI_BSPLINE_HPP
