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


#ifndef CASADI_CONVEXIFY_HPP
#define CASADI_CONVEXIFY_HPP

#include "mx_node.hpp"
#include <map>
#include <stack>

/// \cond INTERNAL

namespace casadi {




  /** \brief Convexify a symmetric matrix
      \author Joris Gillis
      \date 2020
  */
  class CASADI_EXPORT Convexify : public MXNode {
  public:

    /// Constructor
    Convexify(const MX& H, const Dict& opts=Dict());

    /// Destructor
    ~Convexify() override {}

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /** \brief  Evaluate symbolically (MX) */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Generate code for the operation */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief Get required length of iw field */
    size_t sz_iw() const override;

    /** \brief Get required length of w field */
    size_t sz_w() const override;

    /** \brief  Print expression */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Get the operation */
    casadi_int op() const override { return OP_CONVEXIFY;}

    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream& s) const override;

    /** \brief Deserialize without type information */
    static MXNode* deserialize(DeserializingStream& s) { return new Convexify(s); }

    struct ConvexifyData convexify_data_;

    static Sparsity setup(ConvexifyData& d, const Sparsity& H,
                      const Dict& opts=Dict(), bool inplace=true);

    static std::string generate(CodeGenerator& g,
      const ConvexifyData &d,
      const std::string& Hin, const std::string& Hout,
      const std::string& iw, const std::string& w);

    static void deserialize(DeserializingStream& s, const std::string& prefix, ConvexifyData& d);
    static void serialize(SerializingStream& s, const std::string& prefix, const ConvexifyData& d);
  protected:

    /** \brief Deserializing constructor */
    explicit Convexify(DeserializingStream& s);
  };


} // namespace casadi
/// \endcond

#endif // CASADI_CONVEXIFY_HPP
