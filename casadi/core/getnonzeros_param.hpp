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


#ifndef CASADI_GETNONZEROS_PARAM_HPP
#define CASADI_GETNONZEROS_PARAM_HPP

#include "mx_node.hpp"
#include <map>
#include <stack>

/// \cond INTERNAL

namespace casadi {
  /** \brief Get nonzeros of a matrix, parametrically

      \author Joris Gillis
      \date 2019

      \identifier{7v} */
  class CASADI_EXPORT GetNonzerosParam : public MXNode {
  public:

    ///@{
    /// Create functions
    static MX create(const MX& x, const MX& nz);
    static MX create(const MX& x, const MX& inner, const Slice& outer);
    static MX create(const MX& x, const Slice& inner, const MX& outer);
    static MX create(const MX& x, const MX& inner, const MX& outer);
    ///@}

    /// Constructor
    GetNonzerosParam(const Sparsity& sp, const MX& y, const MX& nz);
    GetNonzerosParam(const Sparsity& sp, const MX& y, const MX& nz, const MX& nz_extra);

    /// Destructor
    ~GetNonzerosParam() override {}

    /** \brief  Propagate sparsity forward

        \identifier{7w} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{7x} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief Get the operation

        \identifier{7y} */
    casadi_int op() const override { return OP_GETNONZEROS_PARAM;}

    /** \brief Deserialize without type information

        \identifier{7z} */
    static MXNode* deserialize(DeserializingStream& s);

  protected:
    /** \brief Deserializing constructor

        \identifier{80} */
    explicit GetNonzerosParam(DeserializingStream& s) : MXNode(s) {}
  };


  class CASADI_EXPORT GetNonzerosParamVector : public GetNonzerosParam {
  public:
    /// Constructor
    GetNonzerosParamVector(const MX& x,
                      const MX& nz) : GetNonzerosParam(nz.sparsity(), x, nz) {}

    /// Destructor
    ~GetNonzerosParamVector() override {}

    /** \brief Calculate forward mode directional derivatives

        \identifier{81} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{82} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /** \brief  Evaluate symbolically (MX)

        \identifier{83} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief  Print expression

        \identifier{84} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Generate code for the operation

        \identifier{85} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief Serialize an object without type information

        \identifier{86} */
    void serialize_body(SerializingStream& s) const override;
    /** \brief Serialize type information

        \identifier{87} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserializing constructor

        \identifier{88} */
    explicit GetNonzerosParamVector(DeserializingStream& s);
  };

  // Specialization of the above when nz_ is a nested Slice
  class CASADI_EXPORT GetNonzerosSliceParam : public GetNonzerosParam {
  public:

    /// Constructor
    GetNonzerosSliceParam(const Sparsity& sp, const MX& x, const Slice& inner,
                      const MX& outer) :
                      GetNonzerosParam(sp, x, outer), inner_(inner) {}

    /// Destructor
    ~GetNonzerosSliceParam() override {}

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /** \brief  Evaluate symbolically (MX)

        \identifier{89} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{8a} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{8b} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief  Print expression

        \identifier{8c} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Generate code for the operation

        \identifier{8d} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** Obtain information about node */
    Dict info() const override { return {{"inner", inner_.info()}}; }

    // Data members
    Slice inner_;

    /** \brief Serialize an object without type information

        \identifier{8e} */
    void serialize_body(SerializingStream& s) const override;
    /** \brief Serialize type information

        \identifier{8f} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserializing constructor

        \identifier{8g} */
    explicit GetNonzerosSliceParam(DeserializingStream& s);
  };

  // Specialization of the above when nz_ is a nested Slice
  class CASADI_EXPORT GetNonzerosParamSlice : public GetNonzerosParam {
  public:

    /// Constructor
    GetNonzerosParamSlice(const Sparsity& sp, const MX& x, const MX& inner,
                      const Slice& outer) :
                      GetNonzerosParam(sp, x, inner), outer_(outer) {}

    /// Destructor
    ~GetNonzerosParamSlice() override {}

    /** \brief Get required length of iw field

        \identifier{8h} */
    size_t sz_iw() const override;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /** \brief  Evaluate symbolically (MX)

        \identifier{8i} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{8j} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{8k} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief  Print expression

        \identifier{8l} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Generate code for the operation

        \identifier{8m} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** Obtain information about node */
    Dict info() const override { return {{"outer", outer_.info()}}; }

    // Data members
    Slice outer_;

    /** \brief Serialize an object without type information

        \identifier{8n} */
    void serialize_body(SerializingStream& s) const override;
    /** \brief Serialize type information

        \identifier{8o} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserializing constructor

        \identifier{8p} */
    explicit GetNonzerosParamSlice(DeserializingStream& s);
  };


  // Specialization of the above when nz_ is a nested Slice
  class CASADI_EXPORT GetNonzerosParamParam : public GetNonzerosParam {
  public:

    /// Constructor
    GetNonzerosParamParam(const Sparsity& sp, const MX& x, const MX& inner,
                      const MX& outer) :
                      GetNonzerosParam(sp, x, inner, outer) {}

    /// Destructor
    ~GetNonzerosParamParam() override {}

    /** \brief Get required length of iw field

        \identifier{8q} */
    size_t sz_iw() const override;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /** \brief  Evaluate symbolically (MX)

        \identifier{8r} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{8s} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{8t} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief  Print expression

        \identifier{8u} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Generate code for the operation

        \identifier{8v} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** Obtain information about node */
    Dict info() const override { return {}; }


    /** \brief Serialize type information

        \identifier{8w} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserializing constructor

        \identifier{8x} */
    explicit GetNonzerosParamParam(DeserializingStream& s);
  };

} // namespace casadi
/// \endcond

#endif // CASADI_GETNONZEROS_HPP
