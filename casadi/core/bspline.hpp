/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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

      \identifier{1yd} */
  class CASADI_EXPORT BSplineCommon : public MXNode {
  public:

    /// Constructor
    BSplineCommon(const std::vector<double>& knots,
            const std::vector<casadi_int>& offset,
            const std::vector<casadi_int>& degree,
            casadi_int m,
            const std::vector<casadi_int>& lookup_mode);

    /// Destructor
    ~BSplineCommon() override {}

    static void prepare(casadi_int m, const std::vector<casadi_int>& offset,
      const std::vector<casadi_int>& degree, casadi_int &coeffs_size,
      std::vector<casadi_int>& coeffs_dims, std::vector<casadi_int>& strides);

    static casadi_int get_coeff_size(casadi_int m, const std::vector<casadi_int>& offset,
      const std::vector<casadi_int>& degree);

    template<class M>
    static M derivative_coeff(casadi_int i,
        const std::vector<double>& knots,
        const std::vector<casadi_int>& offset,
        const std::vector<casadi_int>& degree,
        const std::vector<casadi_int>& coeffs_dims,
        const M& coeffs,
        std::vector< std::vector<double> >& new_knots,
        std::vector<casadi_int>& new_degree);

    std::vector<double> knots_;
    std::vector<casadi_int> offset_;
    std::vector<casadi_int> degree_;
    casadi_int m_;
    std::vector<casadi_int> lookup_mode_;

    // Derived fields
    std::vector<casadi_int> strides_;
    std::vector<casadi_int> coeffs_dims_;
    casadi_int coeffs_size_;

    /** \brief Jacobian
     * 
     * Derivatives are computed by transforming the coefficient matrix
     * This is efficient

        \identifier{1ye} */
    mutable MX jac_cache_;

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    /// Mutex for thread safety
    mutable std::mutex jac_cache_mtx_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

    virtual MX jac_cached() const = 0;

    /** \brief Get required length of iw field

        \identifier{1yf} */
    static size_t n_iw(const std::vector<casadi_int> &degree);

    /** \brief Get required length of w field

        \identifier{1yg} */
    static size_t n_w(const std::vector<casadi_int> &degree);

    /** \brief Get required length of iw field

        \identifier{1yh} */
    size_t sz_iw() const override;

    /** \brief Get required length of w field

        \identifier{1yi} */
    size_t sz_w() const override;

    /** \brief Get the operation

        \identifier{1yj} */
    casadi_int op() const override { return OP_BSPLINE;}

    /** \brief Calculate forward mode directional derivatives

        \identifier{1yk} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{1yl} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief Generate code for the operation

        \identifier{1ym} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res,
                  const std::vector<bool>& arg_is_ref,
                  std::vector<bool>& res_is_ref) const override;

    /** \brief Generate code for the operation

        \identifier{1yn} */
    virtual std::string generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<bool>& arg_is_ref) const = 0;

    /** \brief Deserialize without type information

        \identifier{1yo} */
    static MXNode* deserialize(DeserializingStream& s);

    template<class T>
    MX jac(const MX& x, const T& coeffs) const;

    /** \brief Serialize an object without type information

        \identifier{1yp} */
    void serialize_body(SerializingStream& s) const override;

  protected:

    /** \brief Deserializing constructor

        \identifier{1yq} */
    explicit BSplineCommon(DeserializingStream& s);

  };

  /** 
   * 
   * y = bspline(position=symbolic(x),coeffs=numeric);
   * 
   * y in R^m
   * x in R^n
   * 
   */
  class CASADI_EXPORT BSpline : public BSplineCommon {
  public:

    static MX create(const MX& x, const std::vector< std::vector<double> >& knots,
          const std::vector<double>& coeffs,
          const std::vector<casadi_int>& degree,
          casadi_int m,
          const Dict& opts);

    /// Constructor
    BSpline(const MX& x, const std::vector<double>& knots,
            const std::vector<casadi_int>& offset,
            const std::vector<double>& coeffs,
            const std::vector<casadi_int>& degree,
            casadi_int m,
            const std::vector<casadi_int>& lookup_mode);

    /// Destructor
    ~BSpline() override {}

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /** \brief  Evaluate symbolically (MX)

        \identifier{1yr} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Generate code for the operation

        \identifier{1ys} */
    std::string generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<bool>& arg_is_ref) const override;

    /** \brief  Print expression

        \identifier{1yt} */
    std::string disp(const std::vector<std::string>& arg) const override;

    // Numeric coefficients
    std::vector<double> coeffs_;

    MX jac_cached() const override;

    /**
     * 
     * y = bspline(position=numeric(x),coeffs);
     * 
     * x in R^(n x N)
     * y in R^(1 x N)
     * 
     * vec(y) = A coeffs
     * 
     */
    static DM dual(const std::vector<double>& x,
          const std::vector< std::vector<double> >& knots,
          const std::vector<casadi_int>& degree,
          const Dict& opts);
    /** \brief Serialize an object without type information

        \identifier{1yu} */
    void serialize_body(SerializingStream& s) const override;
    /** \brief Serialize type information

        \identifier{1yv} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserializing constructor

        \identifier{1yw} */
    explicit BSpline(DeserializingStream& s);
  };

  // Symbolic coefficients
  class CASADI_EXPORT BSplineParametric : public BSplineCommon {
  public:
    static MX create(const MX& x, const MX& coeffs,
          const std::vector< std::vector<double> >& knots,
          const std::vector<casadi_int>& degree,
          casadi_int m,
          const Dict& opts);

    /// Constructor
    BSplineParametric(const MX& x, const MX& coeffs,
            const std::vector<double>& knots,
            const std::vector<casadi_int>& offset,
            const std::vector<casadi_int>& degree,
            casadi_int m,
            const std::vector<casadi_int>& lookup_mode);

    /// Destructor
    ~BSplineParametric() override {}

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /** \brief  Evaluate symbolically (MX)

        \identifier{1yx} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    MX jac_cached() const override;

    /** \brief Generate code for the operation

        \identifier{1yy} */
    std::string generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<bool>& arg_is_ref) const override;

    /** \brief  Print expression

        \identifier{1yz} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Serialize type information

        \identifier{1z0} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserializing constructor

        \identifier{1z1} */
    explicit BSplineParametric(DeserializingStream& s) : BSplineCommon(s) {}
  };

} // namespace casadi
/// \endcond

#endif // CASADI_BSPLINE_HPP
