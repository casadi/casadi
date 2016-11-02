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


#ifndef CASADI_ROOTFINDER_IMPL_HPP
#define CASADI_ROOTFINDER_IMPL_HPP

#include "rootfinder.hpp"
#include "oracle_function.hpp"
#include "plugin_interface.hpp"

/// \cond INTERNAL
namespace casadi {

  /** \brief Integrator memory */
  struct CASADI_EXPORT RootfinderMemory : public OracleMemory {
    // Inputs
    const double** iarg;

    // Outputs
    double** ires;
  };

  /// Internal class
  class CASADI_EXPORT
  Rootfinder : public OracleFunction, public PluginInterface<Rootfinder> {
  public:
    /** \brief Constructor
     *
     * \param f   Function mapping from (n+1) inputs to 1 output.
     */
    Rootfinder(const std::string& name, const Function& oracle);

    /// Destructor
    virtual ~Rootfinder() = 0;
    ///@{
    /** \brief Number of function inputs and outputs */
    virtual size_t get_n_in() { return oracle_.n_in();}
    virtual size_t get_n_out() { return oracle_.n_out();}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    virtual Sparsity get_sparsity_in(int i) { return oracle_.sparsity_in(i);}
    virtual Sparsity get_sparsity_out(int i) { return oracle_.sparsity_out(i);}
    /// @}

    ///@{
    /** \brief Names of function input and outputs */
    virtual std::string get_name_in(int i) { return oracle_.name_in(i);}
    virtual std::string get_name_out(int i) { return oracle_.name_out(i);}
    /// @}

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    /// Initialize
    virtual void init(const Dict& opts);

    /** \brief Initalize memory block */
    virtual void init_memory(void* mem) const;

    /** \brief Set the (persistent) work vectors */
    virtual void set_work(void* mem, const double**& arg, double**& res,
                          int*& iw, double*& w) const;

    // Evaluate numerically
    virtual void eval(void* mem, const double** arg, double** res, int* iw, double* w) const;

    // Solve the NLP
    virtual void solve(void* mem) const = 0;

    /** \brief  Propagate sparsity forward */
    virtual void sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief  Propagate sparsity backwards */
    virtual void sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    ///@{
    /// Is the class able to propagate seeds through the algorithm?
    virtual bool has_spfwd() const { return true;}
    virtual bool has_sprev() const { return true;}
    ///@}

    ///@{
    /** \brief Generate a function that calculates \a nfwd forward derivatives */
    virtual Function get_forward(const std::string& name, int nfwd,
                                 const std::vector<std::string>& i_names,
                                 const std::vector<std::string>& o_names,
                                 const Dict& opts);
    virtual int get_n_forward() const { return 64;}
    ///@}

    ///@{
    /** \brief Generate a function that calculates \a nadj adjoint derivatives */
    virtual Function get_reverse(const std::string& name, int nadj,
                                 const std::vector<std::string>& i_names,
                                 const std::vector<std::string>& o_names,
                                 const Dict& opts);
    virtual int get_n_reverse() const { return 64;}
    ///@}

    /** \brief Create call to (cached) derivative function, forward mode  */
    virtual void forward(const std::vector<MX>& arg, const std::vector<MX>& res,
                         const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens,
                         bool always_inline, bool never_inline);

    /** \brief Create call to (cached) derivative function, reverse mode  */
    virtual void reverse(const std::vector<MX>& arg, const std::vector<MX>& res,
                         const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens,
                         bool always_inline, bool never_inline);

    /// Number of equations
    int n_;

    /// Linear solver
    Linsol linsol_;
    Sparsity sp_jac_;

    /// Constraints on decision variables
    std::vector<int> u_c_;

    /// Indices of the input and output that correspond to the actual root-finding
    int iin_, iout_;

    // Creator function for internal class
    typedef Rootfinder* (*Creator)(const std::string& name, const Function& oracle);

    // No static functions exposed
    struct Exposed{ };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Short name
    static std::string shortname() { return "rootfinder";}

    /// Infix
    static const std::string infix_;
  };



} // namespace casadi
/// \endcond

#endif // CASADI_ROOTFINDER_IMPL_HPP
