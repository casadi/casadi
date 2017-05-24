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


#ifndef CASADI_NLPSOL_IMPL_HPP
#define CASADI_NLPSOL_IMPL_HPP

#include "nlpsol.hpp"
#include "oracle_function.hpp"
#include "plugin_interface.hpp"


/// \cond INTERNAL
namespace casadi {

  /** \brief Integrator memory */
  struct CASADI_EXPORT NlpsolMemory : public OracleMemory {
    // Inputs
    const double *x0, *p, *lbx, *ubx, *lbg, *ubg, *lam_x0, *lam_g0;

    // Outputs
    double *x, *f, *g, *lam_x, *lam_g, *lam_p;

    // number of iterations
    int n_iter;
  };

  /** \brief NLP solver storage class

      @copydoc Nlpsol_doc
      \author Joel Andersson
      \date 2010-2013
  */
  class CASADI_EXPORT
  Nlpsol : public OracleFunction, public PluginInterface<Nlpsol> {
  public:
    /// Number of variables
    int nx_;

    /// Number of constraints
    int ng_;

    /// Number of parameters
    int np_;

    /// callback function, executed at each iteration
    Function fcallback_;

    /// Execute the callback function only after this amount of iterations
    int callback_step_;

    // Evaluation errors are fatal
    bool eval_errors_fatal_;

    // Warn if initial bounds are violated
    bool warn_initial_bounds_;

    // Ignore errors in the iteration callbacks
    bool iteration_callback_ignore_errors_;

    /// Which variables are discrete?
    std::vector<bool> discrete_;

    // Mixed integer problem?
    bool mi_;

    /// Constructor
    Nlpsol(const std::string& name, const Function& oracle);

    /// Destructor
    ~Nlpsol() override = 0;

    /** \brief Get type name */
    std::string type_name() const override {
      return std::string("nlpsol_") + plugin_name();
    }

    ///@{
    /** \brief Number of function inputs and outputs */
    size_t get_n_in() override { return NLPSOL_NUM_IN;}
    size_t get_n_out() override { return NLPSOL_NUM_OUT;}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    Sparsity get_sparsity_in(int i) override;
    Sparsity get_sparsity_out(int i) override;
    /// @}

    ///@{
    /** \brief Names of function input and outputs */
    std::string get_name_in(int i) override { return nlpsol_in(i);}
    std::string get_name_out(int i) override { return nlpsol_out(i);}
    /// @}

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /// Initialize
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_memory() const override { return new NlpsolMemory();}

    /** \brief Free memory block */
    void free_memory(void *mem) const override { delete static_cast<NlpsolMemory*>(mem);}

    /** \brief Initalize memory block */
    void init_memory(void* mem) const override;

    /** \brief Check if the inputs correspond to a well-posed problem */
    virtual void checkInputs(void* mem) const;

    /** \brief Get default input value */
    double default_in(int ind) const override { return nlpsol_default_in(ind);}

    /// Can discrete variables be treated
    virtual bool integer_support() const { return false;}

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                          int*& iw, double*& w) const override;

    // Evaluate numerically
    void eval(void* mem, const double** arg, double** res, int* iw, double* w) const override;

    // Solve the NLP
    virtual void solve(void* mem) const = 0;

    // Creator function for internal class
    typedef Nlpsol* (*Creator)(const std::string& name, const Function& oracle);

    // No static functions exposed
    struct Exposed{ };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Infix
    static const std::string infix_;

    /// Short name
    static std::string shortname() { return "nlpsol";}

    // Get reduced Hessian
    virtual DM getReducedHessian();

    /// Read options from parameter xml
    virtual void setOptionsFromFile(const std::string & file);

    /// WORKAROUND: Add an element to an std::vector stored in a GenericType:
    template<typename Type> static void append_to_vec(GenericType& t, Type el) {
      std::vector<Type> v = t;
      v.push_back(el);
      t = v;
    }

    /// Convert dictionary to Problem
    template<typename XType>
      static Function map2problem(const std::map<std::string, XType>& d);
  };

  template<typename XType>
  Function Nlpsol::map2problem(const std::map<std::string, XType>& d) {
    std::vector<XType> nl_in(NL_NUM_IN), nl_out(NL_NUM_OUT);
    for (auto&& i : d) {
      if (i.first=="x") {
        nl_in[NL_X]=i.second;
      } else if (i.first=="p") {
        nl_in[NL_P]=i.second;
      } else if (i.first=="f") {
        nl_out[NL_F]=i.second;
      } else if (i.first=="g") {
        nl_out[NL_G]=i.second;
      } else {
        casadi_error("No such field: " + i.first);
      }
    }
    if (nl_out[NL_F].is_empty()) nl_out[NL_F] = 0;
    if (nl_out[NL_G].is_empty()) nl_out[NL_G] = XType(0, 1);
    return Function("nlp", nl_in, nl_out, NL_INPUTS, NL_OUTPUTS);
  }

} // namespace casadi
/// \endcond
#endif // CASADI_NLPSOL_IMPL_HPP
