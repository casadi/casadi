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


#ifndef CASADI_MX_LAYOUT_HPP
#define CASADI_MX_LAYOUT_HPP

#include "mx_node.hpp"

/// \cond INTERNAL
namespace casadi {

  /** \brief MX layout
   * 
   *  y = x.permute_layout(target,range())  <=>
   *          find y such that
   *            interpret(x,x.layout()) == interpret(y,target)
   * 
   * Corollary: if x.layout()==target : identity transform
   * 
   * Lot of examples in unittests
   * 
   * It seems we can get all desired functionality with <perm>. 
   * E.g. we can obtain all tensor-transposes with <perm>=range()
   * 
   * It's not that <perm> is deemed redundant without loss of generality,
   * it's that we currently have no idea what it should mean or when it should be useful
   *
   * \author Joris Gillis
   * \date 2021
  */
  class CASADI_EXPORT PermuteLayout : public MXNode {
  public:
    /// Create function
    static MX create(const MX& x, const Relayout& relay);

    /** \brief  Constructor */
    PermuteLayout(const MX& x, const Relayout& relay);

    /** \brief  Destructor */
    ~PermuteLayout() override {}

    /// Evaluate the function (template)
    //template<typename T>
    //int eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Propagate sparsity forward */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Evaluate symbolically (MX) */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res,
      const std::vector<bool>& unique) const override;

    /** \brief Calculate forward mode directional derivatives */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                          std::vector<std::vector<MX> >& asens) const override;

    /// 1-norm
    MX get_permute_layout(const Relayout& relay) const override;

    /** \brief Generate code for the operation */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res,
                  const std::vector<bool>& arg_is_ref,
                  std::vector<bool>& res_is_ref,
                  bool prefer_inline=false) const override;

    /** \brief  Print expression */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Get the operation */
    casadi_int op() const override { return OP_PERMUTE_LAYOUT;}

    /** \brief Get required length of w field */
    size_t sz_iw() const override;

    /** \brief Deserialize without type information */
    static MXNode* deserialize(DeserializingStream& s) { return new PermuteLayout(s); }

    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream& s) const override;

    MX get_nzref(const Sparsity& sp, const std::vector<casadi_int>& nz,
      bool unique) const override;

    MX as_nzref() const;

  private:
    Relayout relay_;


  protected:
    /** \brief Deserializing constructor */
    explicit PermuteLayout(DeserializingStream& s);
  };

} // namespace casadi

/// \endcond

#endif // CASADI_MX_LAYOUT_HPP
