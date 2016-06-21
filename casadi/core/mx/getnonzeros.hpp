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


#ifndef CASADI_GETNONZEROS_HPP
#define CASADI_GETNONZEROS_HPP

#include "mx_node.hpp"
#include <map>
#include <stack>

/// \cond INTERNAL

namespace casadi {
  /** \brief Get nonzeros of a matrix
      \author Joel Andersson
      \date 2013
  */
  class CASADI_EXPORT GetNonzeros : public MXNode {
  public:

    /// Constructor
    GetNonzeros(const Sparsity& sp, const MX& x);

    /// Destructor
    virtual ~GetNonzeros() {}

    /** \brief  Evaluate symbolically (MX) */
    virtual void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res);

    /** \brief Calculate forward mode directional derivatives */
    virtual void evalFwd(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void evalAdj(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens);

    /// Get an IM representation of a GetNonzeros or SetNonzeros node
    virtual Matrix<int> mapping() const;

    /// Get all the nonzeros
    virtual std::vector<int> all() const = 0;

    /** \brief Get the operation */
    virtual int op() const { return OP_GETNONZEROS;}

    /// Get the nonzeros of matrix
    virtual MX getGetNonzeros(const Sparsity& sp, const std::vector<int>& nz) const;
  };

  class CASADI_EXPORT GetNonzerosVector : public GetNonzeros {
  public:
    /// Constructor
    GetNonzerosVector(const Sparsity& sp, const MX& x,
                      const std::vector<int>& nz) : GetNonzeros(sp, x), nz_(nz) {}

    /// Destructor
    virtual ~GetNonzerosVector() {}

    /// Get all the nonzeros
    virtual std::vector<int> all() const { return nz_;}

    /** \brief  Propagate sparsity forward */
    virtual void sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief  Propagate sparsity backwards */
    virtual void sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /// Evaluate the function (template)
    template<typename T>
    void evalGen(const T* const* arg, T* const* res, int* iw, T* w) const;

    /// Evaluate the function numerically
    virtual void eval(const double** arg, double** res, int* iw, double* w, int mem) const;

    /// Evaluate the function symbolically (SX)
    virtual void eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem);

    /** \brief  Print expression */
    virtual std::string print(const std::vector<std::string>& arg) const;

    /** \brief Generate code for the operation */
    virtual void generate(CodeGenerator& g, const std::string& mem,
                          const std::vector<int>& arg, const std::vector<int>& res) const;

    /** \brief Check if two nodes are equivalent up to a given depth */
    virtual bool is_equal(const MXNode* node, int depth) const;

    /// Operation sequence
    std::vector<int> nz_;
  };

  // Specialization of the above when nz_ is a Slice
  class CASADI_EXPORT GetNonzerosSlice : public GetNonzeros {
  public:

    /// Constructor
    GetNonzerosSlice(const Sparsity& sp, const MX& x, const Slice& s) : GetNonzeros(sp, x), s_(s) {}

    /// Destructor
    virtual ~GetNonzerosSlice() {}

    /// Get all the nonzeros
    virtual std::vector<int> all() const { return s_.all(s_.stop);}

    /// Check if the instance is in fact an identity mapping (that can be simplified)
    bool is_identity() const;

    /// Simplify
    virtual void simplifyMe(MX& ex);

    /** \brief  Propagate sparsity forward */
    virtual void sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief  Propagate sparsity backwards */
    virtual void sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /// Evaluate the function (template)
    template<typename T>
    void evalGen(const T* const* arg, T* const* res, int* iw, T* w) const;

    /// Evaluate the function numerically
    virtual void eval(const double** arg, double** res, int* iw, double* w, int mem) const;

    /// Evaluate the function symbolically (SX)
    virtual void eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem);

    /** \brief  Print expression */
    virtual std::string print(const std::vector<std::string>& arg) const;

    /** \brief Generate code for the operation */
    virtual void generate(CodeGenerator& g, const std::string& mem,
                          const std::vector<int>& arg, const std::vector<int>& res) const;

    /** \brief Check if two nodes are equivalent up to a given depth */
    virtual bool is_equal(const MXNode* node, int depth) const;

    // Data member
    Slice s_;
  };

  // Specialization of the above when nz_ is a nested Slice
  class CASADI_EXPORT GetNonzerosSlice2 : public GetNonzeros {
  public:

    /// Constructor
    GetNonzerosSlice2(const Sparsity& sp, const MX& x, const Slice& inner,
                      const Slice& outer) : GetNonzeros(sp, x), inner_(inner), outer_(outer) {}

    /// Destructor
    virtual ~GetNonzerosSlice2() {}

    /// Get all the nonzeros
    virtual std::vector<int> all() const { return inner_.all(outer_, outer_.stop);}

    /** \brief  Propagate sparsity forward */
    virtual void sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief  Propagate sparsity backwards */
    virtual void sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /// Evaluate the function (template)
    template<typename T>
    void evalGen(const T* const* arg, T* const* res, int* iw, T* w) const;

    /// Evaluate the function numerically
    virtual void eval(const double** arg, double** res, int* iw, double* w, int mem) const;

    /// Evaluate the function symbolically (SX)
    virtual void eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem);

    /** \brief  Print expression */
    virtual std::string print(const std::vector<std::string>& arg) const;

    /** \brief Generate code for the operation */
    virtual void generate(CodeGenerator& g, const std::string& mem,
                          const std::vector<int>& arg, const std::vector<int>& res) const;

    /** \brief Check if two nodes are equivalent up to a given depth */
    virtual bool is_equal(const MXNode* node, int depth) const;

    // Data members
    Slice inner_, outer_;
  };


} // namespace casadi
/// \endcond

#endif // CASADI_GETNONZEROS_HPP
