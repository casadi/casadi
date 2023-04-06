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


#ifndef CASADI_NORM_HPP
#define CASADI_NORM_HPP

#include "mx_node.hpp"

/// \cond INTERNAL
namespace casadi {

  /** \brief Matrix and vector norms

      \author Joel Andersson
      \date 2010-2013

      \identifier{1en} */
  class CASADI_EXPORT Norm : public MXNode {
  public:

    /** \brief  Constructor

        \identifier{1eo} */
    explicit Norm(const MX& x);

    /** \brief  Destructor

        \identifier{1ep} */
    ~Norm() override {}

  protected:
    /** \brief Deserializing constructor

        \identifier{1eq} */
    explicit Norm(DeserializingStream& s) : MXNode(s) {}
  };

  /** \brief Represents a Frobenius norm

      \author Joel Andersson
      \date 2010-2013

      \identifier{1er} */
  class CASADI_EXPORT NormF : public Norm {
  public:

    /** \brief  Constructor

        \identifier{1es} */
    explicit NormF(const MX& x) : Norm(x) {}

    /** \brief  Destructor

        \identifier{1et} */
    ~NormF() override {}

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Evaluate symbolically (MX)

        \identifier{1eu} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{1ev} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{1ew} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief Generate code for the operation

        \identifier{1ex} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief  Print expression

        \identifier{1ey} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Get the operation

        \identifier{1ez} */
    casadi_int op() const override { return OP_NORMF;}

    /** \brief Deserialize without type information

        \identifier{1f0} */
    static MXNode* deserialize(DeserializingStream& s) { return new NormF(s); }

  protected:
    /** \brief Deserializing constructor

        \identifier{1f1} */
    explicit NormF(DeserializingStream& s) : Norm(s) {}
  };

  /** \brief Represents a 2-norm (spectral norm)

      \author Joel Andersson
      \date 2010-2013

      \identifier{1f2} */
  class CASADI_EXPORT Norm2 : public Norm {
  public:

    /** \brief  Constructor

        \identifier{1f3} */
    explicit Norm2(const MX& x): Norm(x) {}

    /** \brief  Destructor

        \identifier{1f4} */
    ~Norm2() override {}

    /** \brief  Print expression

        \identifier{1f5} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Get the operation

        \identifier{1f6} */
    casadi_int op() const override { return OP_NORM2;}

    /** \brief Deserialize without type information

        \identifier{1f7} */
    static MXNode* deserialize(DeserializingStream& s) { return new Norm2(s); }

  protected:
    /** \brief Deserializing constructor

        \identifier{1f8} */
    explicit Norm2(DeserializingStream& s) : Norm(s) {}
  };

  /** \brief 1-norm

      \author Joel Andersson
      \date 2010-2013

      \identifier{1f9} */
  class CASADI_EXPORT Norm1 : public Norm {
  public:

    /** \brief  Constructor

        \identifier{1fa} */
    explicit Norm1(const MX& x) : Norm(x) {}

    /** \brief  Destructor

        \identifier{1fb} */
    ~Norm1() override {}

    /** \brief  Print expression

        \identifier{1fc} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Get the operation

        \identifier{1fd} */
    casadi_int op() const override { return OP_NORM1;}

    /** \brief Deserialize without type information

        \identifier{1fe} */
    static MXNode* deserialize(DeserializingStream& s) { return new Norm1(s); }

  protected:
    /** \brief Deserializing constructor

        \identifier{1ff} */
    explicit Norm1(DeserializingStream& s) : Norm(s) {}
  };

  /** \brief Represents an infinity-norm operation on a MX

      \author Joel Andersson
      \date 2010

      \identifier{1fg} */
  class CASADI_EXPORT NormInf : public Norm {
  public:

    /** \brief  Constructor

        \identifier{1fh} */
    explicit NormInf(const MX& x) : Norm(x) {}

    /** \brief  Destructor

        \identifier{1fi} */
    ~NormInf() override {}

    /** \brief  Print expression

        \identifier{1fj} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Get the operation

        \identifier{1fk} */
    casadi_int op() const override { return OP_NORMINF;}

    /** \brief Deserialize without type information

        \identifier{1fl} */
    static MXNode* deserialize(DeserializingStream& s) { return new NormInf(s); }

  protected:
    /** \brief Deserializing constructor

        \identifier{1fm} */
    explicit NormInf(DeserializingStream& s) : Norm(s) {}
  };

} // namespace casadi

/// \endcond

#endif // CASADI_NORM_HPP
