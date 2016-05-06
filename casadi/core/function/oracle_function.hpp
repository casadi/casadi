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
#include "../timing.hpp"

/// \cond INTERNAL
namespace casadi {

  /** \brief Function memory with temporary work vectors */
  struct CASADI_EXPORT OracleMemory {
    // Work vectors
    const double** arg;
    double** res;
    int* iw;
    double* w;

    // Function specific statistics
    std::map<std::string, FStats> fstats;
  };

  /** \brief Base class for functions that perform calculation with an oracle
      \author Joel Andersson
      \date 2016
  */
  class CASADI_EXPORT OracleFunction : public FunctionInternal {
  protected:
    /// Oracle: Used to generate other functions
    Function oracle_;

    // All NLP functions
    std::map<std::string, Function> all_functions_;
  public:
    /** \brief  Constructor */
    OracleFunction(const std::string& name, const Function& oracle);

    /** \brief  Destructor */
    virtual ~OracleFunction() = 0;

    /** \brief Get oracle */
    virtual const Function& oracle() const { return oracle_;}

    // Replace MX oracle with SX oracle?
    void expand();

    /** Create an oracle function */
    Function
    create_function(const std::string& fname,
                    const std::vector<std::string>& s_in,
                    const std::vector<std::string>& s_out,
                    const Function::AuxOut& aux=Function::AuxOut(),
                    const Dict& opts=Dict());

    /** Register the function for evaluation and statistics gathering */
    void set_function(const std::string& fname, const Function& fcn);

    // Calculate an oracle function
    int calc_function(OracleMemory* m, const std::string& fcn,
                      std::initializer_list<const double*> arg,
                      std::initializer_list<double*> res) const;

    // Get list of dependency functions
    virtual std::vector<std::string> get_function() const;

    // Get a dependency function
    virtual const Function& get_function(const std::string &name) const;

    // Check if a particular dependency exists
    virtual bool has_function(const std::string& fname) const;

    /** \brief Export / Generate C code for the generated functions */
    virtual void generate_dependencies(const std::string& fname, const Dict& opts);

    /** \brief Initalize memory block */
    virtual void init_memory(void* mem) const;

    /// Print statistics
    void print_fstats(const OracleMemory* m) const;

    /// Get all statistics
    virtual Dict get_stats(void* mem) const;
  };

} // namespace casadi

/// \endcond

#endif // CASADI_ORACLE_FUNCTION_HPP
