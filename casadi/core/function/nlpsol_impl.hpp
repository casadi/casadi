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
#include "function_internal.hpp"
#include "plugin_interface.hpp"
#include "../timing.hpp"


/// \cond INTERNAL
namespace casadi {

  /** \brief Integrator memory */
  struct CASADI_EXPORT NlpsolMemory : public WorkMemory {
    // Outputs
    double *x, *f, *g, *lam_x, *lam_g, *lam_p;

    // Inputs
    const double *x0, *p, *lbx, *ubx, *lbg, *ubg, *lam_x0, *lam_g0;

    // Function specific statistics
    std::map<std::string, FStats> fstats;

    // number of iterations
    int n_iter;
  };

  /** \brief NLP solver storage class

      @copydoc Nlpsol_doc
      \author Joel Andersson
      \date 2010-2013
  */
  class CASADI_EXPORT
  Nlpsol : public FunctionInternal, public PluginInterface<Nlpsol> {
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

    // Print timing statistics
    bool print_time_;

    // All NLP functions
    std::vector<Function> all_functions_;

  private:
    /// The NLP
    Oracle* nlp_;

  public:
    /// Constructor
    Nlpsol(const std::string& name, Oracle* nlp);

    /// Destructor
    virtual ~Nlpsol() = 0;

    ///@{
    /** \brief Number of function inputs and outputs */
    virtual size_t get_n_in() { return NLPSOL_NUM_IN;}
    virtual size_t get_n_out() { return NLPSOL_NUM_OUT;}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    virtual Sparsity get_sparsity_in(int i);
    virtual Sparsity get_sparsity_out(int i);
    /// @}

    ///@{
    /** \brief Names of function input and outputs */
    virtual std::string get_name_in(int i) { return nlpsol_in(i);}
    virtual std::string get_name_out(int i) { return nlpsol_out(i);}
    /// @}

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    /// Initialize
    virtual void init(const Dict& opts);

    /** \brief Create memory block */
    virtual void* alloc_memory() const { return new NlpsolMemory();}

    /** \brief Free memory block */
    virtual void free_memory(void *mem) const { delete static_cast<NlpsolMemory*>(mem);}

    /** \brief Initalize memory block */
    virtual void init_memory(void* mem) const;

    /** \brief Check if the inputs correspond to a well-posed problem */
    virtual void checkInputs(void* mem) const;

    /** \brief Get default input value */
    virtual double default_in(int ind) const { return nlpsol_default_in(ind);}

    /** \brief Set the (persistent) work vectors */
    virtual void set_work(void* mem, const double**& arg, double**& res,
                          int*& iw, double*& w) const;

    /** \brief Set the (temporary) work vectors */
    virtual void set_temp(void* mem, const double** arg, double** res,
                          int* iw, double* w) const;

    /** Create an NLP function */
    Function create_function(const std::string& fname,
                             const std::vector<std::string>& s_in,
                             const std::vector<std::string>& s_out,
                             const std::vector<LinComb>& lincomb=std::vector<LinComb>(),
                             const Dict& opts=Dict(), bool reg=true);

    /** Register the function for evaluation and statistics gathering */
    void register_function(const Function& fcn);

    /** \brief Which variables enter nonlinearly */
    std::vector<bool> nl_var(const std::string& s_in,
                             const std::vector<std::string>& s_out) const;

    // Evaluate numerically
    virtual void eval(void* mem, const double** arg, double** res, int* iw, double* w) const;

    // Solve the NLP
    virtual void solve(void* mem) const = 0;

    // Creator function for internal class
    typedef Nlpsol* (*Creator)(const std::string& name, Oracle* nlp);

    // No static functions exposed
    struct Exposed{ };

    // Calculate an oracle function
    int calc_function(NlpsolMemory* m, const Function& fcn,
                      std::initializer_list<const double*> arg,
                      std::initializer_list<double*> res) const;

   /// Print statistics
   void print_fstats(const NlpsolMemory* m) const;

   /// Get all statistics
   virtual Dict get_stats(void* mem) const;

    /** \brief Export / Generate C code for the dependency function */
    virtual void generate_dependencies(const std::string& fname, const Dict& opts);

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
      static Oracle* map2problem(const std::map<std::string, XType>& d);

    /// Get the (legacy) dae forward function
    template<typename XType>
      static Oracle* fun2problem(Function nlp);
  };

  template<typename XType>
  Oracle* Nlpsol::map2problem(const std::map<std::string, XType>& d) {
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
    return Oracle::construct(nl_in, nl_out, NL_INPUTS, NL_OUTPUTS);
  }

  template<typename XType>
  Oracle* Nlpsol::fun2problem(Function nlp) {
    std::vector<XType> nl_in = XType::get_input(nlp);
    std::vector<XType> nl_out = nlp(nl_in);
    return Oracle::construct(nl_in, nl_out, NL_INPUTS, NL_OUTPUTS);
  }

} // namespace casadi
/// \endcond
#endif // CASADI_NLPSOL_IMPL_HPP
