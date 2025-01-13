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


#ifndef CASADI_SX_FUNCTION_HPP
#define CASADI_SX_FUNCTION_HPP

#include "x_function.hpp"

/// \cond INTERNAL

namespace casadi {
  /** \brief  An atomic operation for the SXElem virtual machine

      \identifier{ua} */
  struct ScalarAtomic {
    int op;     /// Operator index
    int i0;
    union {
      double d;
      struct { int i1, i2; };
    };
  };

/** \brief  Internal node class for SXFunction

    Do not use any internal class directly - always use the public Function
    \author Joel Andersson
    \date 2010-2015

    \identifier{ub} */
class CASADI_EXPORT SXFunction :
        public XFunction<SXFunction, Matrix<SXElem>, SXNode>{
  public:
    /** \brief Constructor

        \identifier{uc} */
    SXFunction(const std::string& name,
               const std::vector<Matrix<SXElem> >& inputv,
               const std::vector<Matrix<SXElem> >& outputv,
               const std::vector<std::string>& name_in,
               const std::vector<std::string>& name_out);

  /** \brief  Destructor

      \identifier{ud} */
  ~SXFunction() override;

  /** \brief  Evaluate numerically, work vectors given

      \identifier{ue} */
  int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

  /** \brief  evaluate symbolically while also propagating directional derivatives

      \identifier{uf} */
  int eval_sx(const SXElem** arg, SXElem** res,
              casadi_int* iw, SXElem* w, void* mem,
              bool always_inline, bool never_inline) const override;

  /** \brief Evaluate symbolically, MX type */
  void eval_mx(const MXVector& arg, MXVector& res,
                bool always_inline, bool never_inline) const override;

  /** Inline calls? */
  bool should_inline(bool with_sx, bool always_inline, bool never_inline) const override;

  /** \brief Calculate forward mode directional derivatives

      \identifier{ug} */
  void ad_forward(const std::vector<std::vector<SX> >& fseed,
                            std::vector<std::vector<SX> >& fsens) const;

  /** \brief Calculate reverse mode directional derivatives

      \identifier{uh} */
  void ad_reverse(const std::vector<std::vector<SX> >& aseed,
                            std::vector<std::vector<SX> >& asens) const;

  /** \brief  Check if smooth

      \identifier{ui} */
  bool is_smooth() const;

  /** \brief  Print the algorithm

      \identifier{uj} */
  void disp_more(std::ostream& stream) const override;

  /** \brief Get type name

      \identifier{uk} */
  std::string class_name() const override {return "SXFunction";}

  /** \brief Check if the function is of a particular type

      \identifier{ul} */
  bool is_a(const std::string& type, bool recursive) const override;

  ///@{
  /** \brief Get function input(s) and output(s)

      \identifier{um} */
  const SX sx_in(casadi_int ind) const override;
  const std::vector<SX> sx_in() const override;
  ///@}

  /// Get free variables (SX)
  std::vector<SX> free_sx() const override {
    std::vector<SX> ret(free_vars_.size());
    std::copy(free_vars_.begin(), free_vars_.end(), ret.begin());
    return ret;
  }

  /** \brief Does the function have free variables

      \identifier{un} */
  bool has_free() const override { return !free_vars_.empty();}

  /** \brief Print free variables

      \identifier{uo} */
  std::vector<std::string> get_free() const override {
    std::vector<std::string> ret;
    for (auto&& e : free_vars_) ret.push_back(e.name());
    return ret;
  }

  /** \brief Hessian (forward over adjoint) via source code transformation

      \identifier{up} */
  SX hess(casadi_int iind=0, casadi_int oind=0);

  /** \brief Get the number of atomic operations

      \identifier{uq} */
  casadi_int n_instructions() const override { return algorithm_.size();}

  /** \brief Get an atomic operation operator index

      \identifier{ur} */
  casadi_int instruction_id(casadi_int k) const override { return algorithm_.at(k).op;}

  /** \brief Get the (integer) input arguments of an atomic operation

      \identifier{us} */
  std::vector<casadi_int> instruction_input(casadi_int k) const override {
    auto e = algorithm_.at(k);
    if (e.op==OP_CALL) {
      const ExtendedAlgEl& m = call_.el[e.i1];
      return vector_static_cast<casadi_int>(m.dep);
    } else if (casadi_math<double>::ndeps(e.op)==2 || e.op==OP_INPUT) {
      return {e.i1, e.i2};
    } else if (casadi_math<double>::ndeps(e.op)==1) {
      return {e.i1};
    } else {
      return {};
    }
  }

  /** \brief Get the floating point output argument of an atomic operation

      \identifier{ut} */
  double instruction_constant(casadi_int k) const override {
    return algorithm_.at(k).d;
  }

  /** \brief Get the (integer) output argument of an atomic operation

      \identifier{uu} */
  std::vector<casadi_int> instruction_output(casadi_int k) const override {
    auto e = algorithm_.at(k);
    if (e.op==OP_CALL) {
      const ExtendedAlgEl& m = call_.el[e.i1];
      return vector_static_cast<casadi_int>(m.res);
    } else if (e.op==OP_OUTPUT) {
      return {e.i0, e.i2};
    } else {
      return {e.i0};
    }
  }

  /** \brief Number of nodes in the algorithm

      \identifier{uv} */
  casadi_int n_nodes() const override { return algorithm_.size() - nnz_out();}

  /** \brief  DATA MEMBERS

      \identifier{uw} */

  /** \brief  An element of the algorithm, namely a binary operation

      \identifier{ux} */
  typedef ScalarAtomic AlgEl;

  /** \brief  An element of the tape

      \identifier{uy} */
  template<typename T>
  struct TapeEl {
    T d[2];
  };

  /** \brief  all binary nodes of the tree in the order of execution

      \identifier{uz} */
  std::vector<AlgEl> algorithm_;

  // Work vector size
  size_t worksize_;

  /// Free variables
  std::vector<SXElem> free_vars_;

  /// The expressions corresponding to each binary operation
  std::vector<SXElem> operations_;

  /// The expressions corresponding to each constant
  std::vector<SXElem> constants_;

  /// Default input values
  std::vector<double> default_in_;

  /// Copy elision per algel
  std::vector<bool> copy_elision_;

    /** \brief Serialize an object without type information

        \identifier{v0} */
  void serialize_body(SerializingStream &s) const override;

  // call node information that won't fit into AlgEl
  struct ExtendedAlgEl {
    ExtendedAlgEl(const Function& fun);
    Function f;
    // Work vector indices of the arguments (cfr AlgEl::arg)
    std::vector<int> dep;
    // Work vector indices of the results (cfr AlgEl::res)
    std::vector<int> res;

    std::vector<int> copy_elision_arg;
    std::vector<int> copy_elision_offset;

    // Following fields are redundant but will increase eval speed
    casadi_int n_dep;
    casadi_int n_res;
    casadi_int f_n_in;
    casadi_int f_n_out;
    std::vector<int> f_nnz_in;
    std::vector<int> f_nnz_out;
  };

  /// Metadata for call nodes
  struct CallInfo {
    // Maximum memory requirements across all call nodes
    size_t sz_arg = 0, sz_res = 0, sz_iw = 0, sz_w = 0;
    size_t sz_w_arg = 0, sz_w_res = 0;
    std::vector<ExtendedAlgEl> el;
  } call_;

    /** \brief Deserialize without type information

        \identifier{v1} */
  static ProtoFunction* deserialize(DeserializingStream& s);

  static std::vector<SX> order(const std::vector<SX>& expr);

  ///@{
  /** \brief Options

      \identifier{v2} */
  static const Options options_;
  const Options& get_options() const override { return options_;}
  ///@}

  /// Reconstruct options dict
  Dict generate_options(const std::string& target="clone") const override;

  /** \brief  Initialize

      \identifier{v3} */
  void init(const Dict& opts) override;

  /** \brief Part of initialize responsible of prepaprign copy elision

      \identifier{29h} */
  void init_copy_elision();

  /** \brief  Get the size of the work vector, for codegen

      \identifier{290} */
  size_t codegen_sz_w(const CodeGenerator& g) const override;

  /** \brief Generate code for the declarations of the C function

      \identifier{v4} */
  void codegen_declarations(CodeGenerator& g) const override;

  /** \brief Generate code for the body of the C function

      \identifier{v5} */
  void codegen_body(CodeGenerator& g) const override;

  /** \brief  Propagate sparsity forward

      \identifier{v6} */
  int sp_forward(const bvec_t** arg, bvec_t** res,
                  casadi_int* iw, bvec_t* w, void* mem) const override;

  /** \brief  Propagate sparsity backwards

      \identifier{v7} */
  int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const override;

  /** *\brief get SX expression associated with instructions

       \identifier{v8} */
  SX instructions_sx() const override;

  // Get all embedded functions, recursively
  void find(std::map<FunctionInternal*, Function>& all_fun, casadi_int max_depth) const override;

  /** \brief Get default input value

      \identifier{v9} */
  double get_default_in(casadi_int ind) const override { return default_in_.at(ind);}

  /** \brief Export function in a specific language

      \identifier{va} */
  void export_code_body(const std::string& lang,
    std::ostream &stream, const Dict& options) const override;

  /// With just-in-time compilation using OpenCL
  bool just_in_time_opencl_;

  /// With just-in-time compilation for the sparsity propagation
  bool just_in_time_sparsity_;

  /// Live variables?
  bool live_variables_;

protected:
  template<typename T>
  void call_fwd(const AlgEl& e, const T** arg, T** res, casadi_int* iw, T* w) const;

  template<typename T>
  void call_rev(const AlgEl& e, T** arg, T** res, casadi_int* iw, T* w) const;

  template<typename T, typename CT>
  void call_setup(const ExtendedAlgEl& m,
    CT*** call_arg, T*** call_res, casadi_int** call_iw, T** call_w, T** nz_in, T** nz_out) const;

  /** \brief Deserializing constructor

      \identifier{vb} */
  explicit SXFunction(DeserializingStream& s);
};


} // namespace casadi

/// \endcond
#endif // CASADI_SX_FUNCTION_HPP
