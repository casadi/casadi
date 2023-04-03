/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
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


#ifndef CASADI_SPLIT_HPP
#define CASADI_SPLIT_HPP

#include "multiple_output.hpp"
#include <map>
#include <stack>

/// \cond INTERNAL

namespace casadi {

  /** \brief Split: Split into multiple expressions splitting the nonzeros

      \author Joel Andersson
      \date 2014

      \identifier{12d} */
  class CASADI_EXPORT Split : public MultipleOutput {
  public:
    /// Constructor
    Split(const MX& x, const std::vector<casadi_int>& offset);

    /// Destructor
    ~Split() override = 0;

    /** \brief  Number of outputs

        \identifier{12e} */
    casadi_int nout() const override { return output_sparsity_.size(); }

    /** \brief  Get the sparsity of output oind

        \identifier{12f} */
    const Sparsity& sparsity(casadi_int oind) const override { return output_sparsity_.at(oind);}

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Propagate sparsity forward

        \identifier{12g} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{12h} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief Generate code for the operation

        \identifier{12i} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** Obtain information about node */
    Dict info() const override;

    // Sparsity pattern of the outputs
    std::vector<casadi_int> offset_;
    std::vector<Sparsity> output_sparsity_;

    /** \brief Serialize an object without type information

        \identifier{12j} */
    void serialize_body(SerializingStream& s) const override;

  protected:
    /** \brief Deserializing constructor

        \identifier{12k} */
    explicit Split(DeserializingStream& s);
  };

  /** \brief Horizontal split, x -> x0, x1, ...

      \author Joel Andersson
      \date 2013

      \identifier{12l} */
  class CASADI_EXPORT Horzsplit : public Split {
  public:

    /// Constructor
    Horzsplit(const MX& x, const std::vector<casadi_int>& offset);

    /// Destructor
    ~Horzsplit() override {}

    /** \brief  Evaluate symbolically (MX)

        \identifier{12m} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{12n} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{12o} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief  Print expression

        \identifier{12p} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Get the operation

        \identifier{12q} */
    casadi_int op() const override { return OP_HORZSPLIT;}

    /// Create a horizontal concatenation node
    MX get_horzcat(const std::vector<MX>& x) const override;

    /** \brief Deserialize without type information

        \identifier{12r} */
    static MXNode* deserialize(DeserializingStream& s) { return new Horzsplit(s); }

  protected:
    /** \brief Deserializing constructor

        \identifier{12s} */
    explicit Horzsplit(DeserializingStream& s) : Split(s) {}
  };

  /** \brief Diag split, x -> x0, x1, ...

      \author Joris Gillis
      \date 2014

      \identifier{12t} */
  class CASADI_EXPORT Diagsplit : public Split {
  public:

    /// Constructor
    Diagsplit(const MX& x,
      const std::vector<casadi_int>& offset1, const std::vector<casadi_int>& offset2);

    /// Destructor
    ~Diagsplit() override {}

    /** \brief  Evaluate symbolically (MX)

        \identifier{12u} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{12v} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{12w} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief  Print expression

        \identifier{12x} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Get the operation

        \identifier{12y} */
    casadi_int op() const override { return OP_DIAGSPLIT;}

    /// Create a diagonal concatenation node
    MX get_diagcat(const std::vector<MX>& x) const override;

    /** \brief Deserialize without type information

        \identifier{12z} */
    static MXNode* deserialize(DeserializingStream& s) { return new Diagsplit(s); }

  protected:
    /** \brief Deserializing constructor

        \identifier{130} */
    explicit Diagsplit(DeserializingStream& s) : Split(s) {}
  };

  /** \brief Vertical split of vectors, x -> x0, x1, ...

      \author Joel Andersson
      \date 2014

      \identifier{131} */
  class CASADI_EXPORT Vertsplit : public Split {
  public:

    /// Constructor
    Vertsplit(const MX& x, const std::vector<casadi_int>& offset);

    /// Destructor
    ~Vertsplit() override {}

    /** \brief  Evaluate symbolically (MX)

        \identifier{132} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{133} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{134} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief  Print expression

        \identifier{135} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Get the operation

        \identifier{136} */
    casadi_int op() const override { return OP_VERTSPLIT;}

    /// Create a vertical concatenation node (vectors only)
    MX get_vertcat(const std::vector<MX>& x) const override;

    /** \brief Deserialize without type information

        \identifier{137} */
    static MXNode* deserialize(DeserializingStream& s) { return new Vertsplit(s); }

  protected:
    /** \brief Deserializing constructor

        \identifier{138} */
    explicit Vertsplit(DeserializingStream& s) : Split(s) {}
  };

} // namespace casadi

/// \endcond

#endif // CASADI_SPLIT_HPP
