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


#ifndef CASADI_FMU_FUNCTION_HPP
#define CASADI_FMU_FUNCTION_HPP

#include "function_internal.hpp"
#include "fmu.hpp"
#include "finite_differences.hpp"

/// \cond INTERNAL

namespace casadi {

// Forward declarations
class DaeBuilderInternal;
class FmuFunction;
struct InputStruct;

// Memory object
struct CASADI_EXPORT FmuMemory : public FunctionMemory {
  // Function object
  const FmuFunction& self;
  // Evaluation inputs
  const double** arg;
  // Evaluation outputs
  double** res;
  // Work vector for star coloring
  casadi_int* star_iw;
  // Extended Jacobian
  double *jac_nz;
  // Extended Hessian
  double *hess_nz;
  // Adjoint seeds, sensitivities being calculated
  double *aseed, *asens, *pert_asens;
  // Memory for Jacobian calculation
  casadi_jac_data<double> d;
  // Instance memory
  void* instance;
  // Additional (slave) memory objects
  std::vector<FmuMemory*> slaves;
  // Input and output buffers
  std::vector<double> ibuf_, obuf_;
  // Seeds, sensitivities
  std::vector<double> seed_, sens_;
  // Which entries have been changed
  std::vector<bool> changed_;
  // Which entries are being requested
  std::vector<bool> requested_;
  // Derivative with respect to
  std::vector<size_t> wrt_;
  // Current known/unknown variables
  std::vector<size_t> id_in_, id_out_;
  // Which perturbations are permitted
  std::vector<bool> in_bounds_;
  // Flip sign?
  std::vector<bool> flip_;
  // Value references
  std::vector<unsigned int> vr_in_, vr_out_;
  // Work vector (reals)
  std::vector<double> v_in_, v_out_, d_in_, d_out_, fd_out_, v_pert_;
  // Constructor
  explicit FmuMemory(const FmuFunction& self) : self(self), instance(nullptr) {}
};

/// Type of parallelization
enum class Parallelization {SERIAL, OPENMP, THREAD, NUMEL};

/// Convert to string
CASADI_EXPORT std::string to_string(Parallelization v);

// Types of inputs
enum class InputType {REG, FWD, ADJ, OUT, ADJ_OUT};

// Input structure
struct CASADI_EXPORT InputStruct {
  // Type of input
  InputType type;
  // Corresponding index in Fmu
  size_t ind;
  // Parse an input string
  static InputStruct parse(const std::string& n, const Fmu* fmu,
    std::vector<std::string>* name_in = 0,
    std::vector<std::string>* name_out = 0);
};

// Types of inputs
enum class OutputType {REG, FWD, ADJ, JAC, JAC_TRANS, JAC_ADJ_OUT, JAC_REG_ADJ, HESS};

// Output structure
struct CASADI_EXPORT OutputStruct {
  // Type of input
  OutputType type;
  // Output index in Fmu
  size_t ind;
  // With-respect-to index in Fmu
  size_t wrt;
  // Selection
  size_t rbegin, rend, cbegin, cend;
  // Parse an output string
  static OutputStruct parse(const std::string& n, const Fmu* fmu,
    std::vector<std::string>* name_in = 0,
    std::vector<std::string>* name_out = 0);
  // Constructor
  OutputStruct() : ind(-1), wrt(-1), rbegin(-1), rend(-1), cbegin(-1), cend(-1) {}
};

// Helper function
CASADI_EXPORT bool has_prefix(const std::string& s);

// Split prefix
CASADI_EXPORT std::string pop_prefix(const std::string& s, std::string* rem = 0);

class CASADI_EXPORT FmuFunction : public FunctionInternal {
 public:
  // FMU (shared between derivative expressions
  Fmu fmu_;

  // Information about function inputs
  std::vector<InputStruct> in_;

  // Information about function outputs
  std::vector<OutputStruct> out_;

  // All Jacobian inputs and outputs
  std::vector<size_t> jac_in_, jac_out_;

  // Nominal values for Jacobian inputs
  std::vector<double> jac_nom_in_;

  // Sparsity of transpose (if needed)
  std::vector<Sparsity> sp_trans_;
  std::vector<casadi_int> sp_trans_map_;

  // What blocks exist?
  bool has_jac_, has_fwd_, has_adj_, has_hess_;

  // Merge with enable_forward_ in FunctionInternal
  bool enable_ad_;

  // Validate derivative calculations: Move to base class?
  bool validate_forward_, validate_hessian_;

  // User-set options
  bool make_symmetric_;
  double step_, abstol_, reltol_;
  bool print_progress_, new_jacobian_, new_forward_, new_hessian_, hessian_coloring_;
  std::string validate_ad_file_;

  // FD method as an enum
  FdMode fd_;

  // Types of parallelization
  Parallelization parallelization_;

  // Stats from initialization
  Dict init_stats_;

  /** \brief Constructor

      \identifier{yh} */
  FmuFunction(const std::string& name, const Fmu& fmu,
      const std::vector<std::string>& name_in,
      const std::vector<std::string>& name_out);

  /** \brief Destructor

      \identifier{yi} */
  ~FmuFunction() override;

  /** \brief Get type name

      \identifier{yj} */
  std::string class_name() const override { return "FmuFunction";}

  ///@{
  /** \brief Options

      \identifier{yk} */
  static const Options options_;
  const Options& get_options() const override { return options_;}
  ///@}

  /// Initialize
  void init(const Dict& opts) override;

  // Identify input and output schemes from FmuFunction inputs and outputs
  static void identify_io(
    std::vector<std::string>* scheme_in,
    std::vector<std::string>* scheme_out,
    const std::vector<std::string>& name_in,
    const std::vector<std::string>& name_out);

  // Get sparsity pattern for extended Jacobian, Hessian
  Sparsity jac_sp_, hess_sp_;

  // Graph coloring
  Sparsity jac_colors_, hess_colors_;

  // Nonlinearly entering variables
  std::vector<casadi_int> nonlin_;

  // Jacobian memory
  casadi_jac_prob<double> p_;

  // Number of parallel tasks
  casadi_int max_jac_tasks_, max_hess_tasks_, max_n_tasks_;

  ///@{
  /** \brief Number of function inputs and outputs

      \identifier{yl} */
  size_t get_n_in() override { return in_.size();}
  size_t get_n_out() override {return out_.size();}
  ///@}

  /// @{
  /** \brief Retreive sparsities

      \identifier{ym} */
  Sparsity get_sparsity_in(casadi_int i) override;
  Sparsity get_sparsity_out(casadi_int i) override;
  /// @}

  /// @{
  /** \brief Retreive nominal values

      \identifier{yn} */
  std::vector<double> get_nominal_in(casadi_int i) const override;
  std::vector<double> get_nominal_out(casadi_int i) const override;
  /// @}

  // Evaluate numerically
  int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

  // Evaluate all tasks numerically, serially or in parallel
  int eval_all(FmuMemory* m, casadi_int n_task,
    bool need_nondiff, bool need_jac, bool need_fwd, bool need_adj, bool need_hess) const;

  // Evaluate numerically, single thread
  int eval_task(FmuMemory* m, casadi_int task, casadi_int n_task,
    bool need_nondiff, bool need_jac, bool need_fwd, bool need_adj, bool need_hess) const;

  // Remove NaNs from Hessian (necessary for star coloring approach)
  void remove_nans(double *hess_nz, casadi_int* iw) const;

  // Check extended Hessian
  void check_hessian(FmuMemory* m, const double *hess_nz, casadi_int* iw) const;

  // Make extended Hessian symmetric
  void make_symmetric(double *hess_nz, casadi_int* iw) const;

  ///@{
  /** \brief Return sparsity of Jacobian of an output respect to an input

      \identifier{yo} */
  bool has_jac_sparsity(casadi_int oind, casadi_int iind) const override;
  Sparsity get_jac_sparsity(casadi_int oind, casadi_int iind, bool symmetric) const override;
  ///@}

  // Are all inputs/outputs regular?
  bool all_regular() const;

  // Are all inputs/outputs vectors (i.e. not Jacobian or Hessian blocks)?
  bool all_vectors() const;

  // Factory
  Function factory(const std::string& name,
    const std::vector<std::string>& s_in,
    const std::vector<std::string>& s_out,
    const Function::AuxOut& aux,
    const Dict& opts) const override;

  ///@{
  /** \brief Full Jacobian

      \identifier{yp} */
  bool has_jacobian() const override;
  Function get_jacobian(const std::string& name,
    const std::vector<std::string>& inames,
    const std::vector<std::string>& onames,
    const Dict& opts) const override;
  ///@}

  ///@{
  /** \brief Return function that calculates forward derivatives

      \identifier{27h} */
  bool has_forward(casadi_int nfwd) const override;
  Function get_forward(casadi_int nfwd, const std::string& name,
    const std::vector<std::string>& inames,
    const std::vector<std::string>& onames,
    const Dict& opts) const override;
  ///@}

  ///@{
  /** \brief Reverse mode AD

      \identifier{yq} */
  bool has_reverse(casadi_int nadj) const override;
  Function get_reverse(casadi_int nadj, const std::string& name,
    const std::vector<std::string>& inames,
    const std::vector<std::string>& onames,
    const Dict& opts) const override;
  ///@}

  /** \brief Create memory block

      \identifier{yr} */
  void* alloc_mem() const override;

  /** \brief Initalize memory block

      \identifier{ys} */
  int init_mem(void* mem) const override;

  /** \brief Free memory block

      \identifier{yt} */
  void free_mem(void *mem) const override;

  /// Get all statistics
  Dict get_stats(void* mem) const override;

  /** \brief Serialize an object without type information

      \identifier{279} */
  void serialize_body(SerializingStream &s) const override;

  /** \brief Deserialize without type information

      \identifier{27a} */
  static ProtoFunction* deserialize(DeserializingStream& s) { return new FmuFunction(s); }

  /** \brief Change option after object creation for debugging

      \identifier{27f} */
  void change_option(const std::string& option_name, const GenericType& option_value) override;

  protected:
    /** \brief Deserializing constructor

        \identifier{27b} */
    explicit FmuFunction(DeserializingStream& s);
};

} // namespace casadi
/// \endcond

#endif // CASADI_FMU_FUNCTION_HPP
