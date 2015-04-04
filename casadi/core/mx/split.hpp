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
  */
  class CASADI_EXPORT Split : public MultipleOutput {
  public:
    /// Constructor
    Split(const MX& x, const std::vector<int>& offset);

    /// Destructor
    virtual ~Split() = 0;

    /** \brief  Number of outputs */
    virtual int nout() const { return output_sparsity_.size(); }

    /** \brief  Get the sparsity of output oind */
    virtual const Sparsity& sparsity(int oind) const { return output_sparsity_.at(oind);}

    /// Evaluate the function (template)
    template<typename T>
    void evalGen(const T* const* arg, T* const* res, int* itmp, T* rtmp);

    /// Evaluate the function numerically
    virtual void evalD(cp_double* input, p_double* output,
                       int* itmp, double* rtmp);

    /// Evaluate the function symbolically (SX)
    virtual void evalSX(cp_SXElement* input, p_SXElement* output,
                            int* itmp, SXElement* rtmp);

    /** \brief  Propagate sparsity forward */
    virtual void spFwd(cp_bvec_t* arg,
                       p_bvec_t* res, int* itmp, bvec_t* rtmp);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdj(p_bvec_t* arg,
                       p_bvec_t* res, int* itmp, bvec_t* rtmp);

    /** \brief Generate code for the operation */
    virtual void generate(std::ostream &stream, const std::vector<int>& arg,
                                   const std::vector<int>& res, CodeGenerator& gen) const;

    // Sparsity pattern of the outputs
    std::vector<int> offset_;
    std::vector<Sparsity> output_sparsity_;
  };

  /** \brief Horizontal split, x -> x0, x1, ...
      \author Joel Andersson
      \date 2013
  */
  class CASADI_EXPORT Horzsplit : public Split {
  public:

    /// Constructor
    Horzsplit(const MX& x, const std::vector<int>& offset);

    /// Destructor
    virtual ~Horzsplit() {}

    /// Clone function
    virtual Horzsplit* clone() const { return new Horzsplit(*this);}

    /// Evaluate the function symbolically (MX)
    virtual void eval(const cpv_MX& input, const pv_MX& output);

    /** \brief Calculate forward mode directional derivatives */
    virtual void evalFwd(const std::vector<cpv_MX>& fwdSeed, const std::vector<pv_MX>& fwdSens);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void evalAdj(const std::vector<pv_MX>& adjSeed, const std::vector<pv_MX>& adjSens);

    /// Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;

    /** \brief Get the operation */
    virtual int getOp() const { return OP_HORZSPLIT;}

    /// Create a horizontal concatenation node
    virtual MX getHorzcat(const std::vector<MX>& x) const;
  };

  /** \brief Diag split, x -> x0, x1, ...
      \author Joris Gillis
      \date 2014
  */
  class CASADI_EXPORT Diagsplit : public Split {
  public:

    /// Constructor
    Diagsplit(const MX& x, const std::vector<int>& offset1, const std::vector<int>& offset2);

    /// Destructor
    virtual ~Diagsplit() {}

    /// Clone function
    virtual Diagsplit* clone() const { return new Diagsplit(*this);}

    /// Evaluate the function symbolically (MX)
    virtual void eval(const cpv_MX& input, const pv_MX& output);

    /** \brief Calculate forward mode directional derivatives */
    virtual void evalFwd(const std::vector<cpv_MX>& fwdSeed, const std::vector<pv_MX>& fwdSens);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void evalAdj(const std::vector<pv_MX>& adjSeed, const std::vector<pv_MX>& adjSens);

    /// Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;

    /** \brief Get the operation */
    virtual int getOp() const { return OP_DIAGSPLIT;}

    /// Create a diagonal concatenation node
    virtual MX getDiagcat(const std::vector<MX>& x) const;
  };

  /** \brief Vertical split of vectors, x -> x0, x1, ...
      \author Joel Andersson
      \date 2014
  */
  class CASADI_EXPORT Vertsplit : public Split {
  public:

    /// Constructor
    Vertsplit(const MX& x, const std::vector<int>& offset);

    /// Destructor
    virtual ~Vertsplit() {}

    /// Clone function
    virtual Vertsplit* clone() const { return new Vertsplit(*this);}

    /// Evaluate the function symbolically (MX)
    virtual void eval(const cpv_MX& input, const pv_MX& output);

    /** \brief Calculate forward mode directional derivatives */
    virtual void evalFwd(const std::vector<cpv_MX>& fwdSeed, const std::vector<pv_MX>& fwdSens);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void evalAdj(const std::vector<pv_MX>& adjSeed, const std::vector<pv_MX>& adjSens);

    /// Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;

    /** \brief Get the operation */
    virtual int getOp() const { return OP_VERTSPLIT;}

    /// Create a vertical concatenation node (vectors only)
    virtual MX getVertcat(const std::vector<MX>& x) const;
  };

} // namespace casadi

/// \endcond

#endif // CASADI_SPLIT_HPP
