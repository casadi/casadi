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


#ifndef CASADI_ORACLE_FUNCTION_HPP
#define CASADI_ORACLE_FUNCTION_HPP

#include "function_internal.hpp"

/// \cond INTERNAL
namespace casadi {

  /** \brief Function memory with temporary work vectors

      \identifier{b} */
  struct CASADI_EXPORT LocalOracleMemory : public FunctionMemory {
    // Work vectors
    const double** arg;
    double** res;
    casadi_int* iw;
    double* w;
  };

  /** \brief Function memory

      \identifier{c} */
  struct CASADI_EXPORT OracleMemory : public FunctionMemory {
    // Work vector aliases for non-threaded convenience
    const double** arg;
    double** res;
    casadi_int* iw;
    double* w;

    std::vector<LocalOracleMemory*> thread_local_mem;
  };

  /** \brief Base class for functions that perform calculation with an oracle

      \author Joel Andersson
      \date 2016

      \identifier{d} */
  class CASADI_EXPORT OracleFunction : public FunctionInternal {
  protected:
    /// Oracle: Used to generate other functions
    Function oracle_;

    /// Options for creating functions
    Dict common_options_;
    Dict specific_options_;

    /// Show evaluation warnings
    bool show_eval_warnings_;

    // Maximum number of threads
    int max_num_threads_;

    // Information about one function
    struct RegFun {
      Function f;
      bool jit;
      Function f_original; // Relevant for jit
      bool monitored = false;
    };

    // All NLP functions
    std::map<std::string, RegFun> all_functions_;

    // Active monitors
    std::vector<std::string> monitor_;

    // Memory stride in case of multipel threads
    size_t stride_arg_, stride_res_, stride_iw_, stride_w_;

  public:
    /** \brief  Constructor

        \identifier{e} */
    OracleFunction(const std::string& name, const Function& oracle);

    /** \brief  Destructor

        \identifier{f} */
    ~OracleFunction() override = 0;

    ///@{
    /** \brief Options

        \identifier{g} */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** Initialize  */
    void init(const Dict& opts) override;

    /// Finalize initialization
    void finalize() override;

    /// Combine results from different threads
    void join_results(OracleMemory* m) const;

    /** \brief Get oracle

        \identifier{h} */
    const Function& oracle() const override { return oracle_;}

    // Replace MX oracle with SX oracle?
    void expand();

    /** Create an oracle function, using a different oracle function
     * Temporary addition to allow transition from one oracle formulation to another.
     */
    Function create_function(const Function& oracle, const std::string& fname,
      const std::vector<std::string>& s_in,
      const std::vector<std::string>& s_out,
      const Function::AuxOut& aux=Function::AuxOut());

    /** Create an oracle function */
    Function create_function(const std::string& fname,
      const std::vector<std::string>& s_in,
      const std::vector<std::string>& s_out,
      const Function::AuxOut& aux=Function::AuxOut());

    /** Create an oracle function as a forward derivative of a different function */
    Function create_forward(const std::string& fname, casadi_int nfwd);

    /** Register the function for evaluation and statistics gathering */
    void set_function(const Function& fcn, const std::string& fname, bool jit=false);

    /** Register the function for evaluation and statistics gathering */
    void set_function(const Function& fcn) { set_function(fcn, fcn.name()); }

    // Calculate an oracle function
    int calc_function(OracleMemory* m, const std::string& fcn,
      const double* const* arg=nullptr, int thread_id=0) const;

    // Forward sparsity propagation through a function
    int calc_sp_forward(const std::string& fcn, const bvec_t** arg, bvec_t** res,
      casadi_int* iw, bvec_t* w) const;

    // Reverse sparsity propagation through a function
    int calc_sp_reverse(const std::string& fcn, bvec_t** arg, bvec_t** res,
      casadi_int* iw, bvec_t* w) const;

    /** \brief Get list of dependency functions

     * -1 Indicates irregularity

        \identifier{i} */
    std::vector<std::string> get_function() const override;

    // Get a dependency function
    const Function& get_function(const std::string &name) const override;

    // Is a function monitored?
    virtual bool monitored(const std::string &name) const;

    // Check if a particular dependency exists
    bool has_function(const std::string& fname) const override;

    /** \brief Export / Generate C code for the generated functions

        \identifier{j} */
    std::string generate_dependencies(const std::string& fname, const Dict& opts) const override;

    /** \brief JIT for dependencies

        \identifier{k} */
    void jit_dependencies(const std::string& fname) override;

    /** \brief Create memory block

        \identifier{l} */
    void* alloc_mem() const override { return new OracleMemory();}

    /** \brief Initalize memory block

        \identifier{m} */
    int local_init_mem(void* mem) const;

    /** \brief Initalize memory block

        \identifier{n} */
    int init_mem(void* mem) const override;

    /** \brief Free memory block

        \identifier{o} */
    void local_free_mem(void *mem) const { delete static_cast<LocalOracleMemory*>(mem);}

   /** \brief Free memory block

       \identifier{p} */
    void free_mem(void *mem) const override;

    /** \brief Set the work vectors

        \identifier{q} */
    void set_temp(void* mem, const double** arg, double** res,
                          casadi_int* iw, double* w) const override;

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    /** \brief Serialize an object without type information

        \identifier{r} */
    void serialize_body(SerializingStream &s) const override;

  protected:
    /** \brief Deserializing constructor

        \identifier{s} */
    explicit OracleFunction(DeserializingStream& s);

  };

} // namespace casadi

/// \endcond

#endif // CASADI_ORACLE_FUNCTION_HPP
