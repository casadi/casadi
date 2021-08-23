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


#ifndef CASADI_SIMULATOR_IMPL_HPP
#define CASADI_SIMULATOR_IMPL_HPP

#include "simulator.hpp"
#include "oracle_function.hpp"
#include "plugin_interface.hpp"

/// \cond INTERNAL

namespace casadi {

  /** \brief Simulator memory */
  struct CASADI_EXPORT SimulatorMemory : public OracleMemory {
  };

  /** \brief Internal storage for simulator related data

      @copydoc DAE_doc
      \author Joel Andersson
      \date 2010
  */
  class CASADI_EXPORT
  Simulator : public OracleFunction, public PluginInterface<Simulator> {
  public:
    /** \brief  Constructor */
    Simulator(const std::string& name, const Function& oracle);

    /** \brief  Destructor */
    ~Simulator() override=0;

    ///@{
    /** \brief Number of function inputs and outputs */
    size_t get_n_in() override { return SIMULATOR_NUM_IN;}
    size_t get_n_out() override { return SIMULATOR_NUM_OUT;}
    ///@}

   /// @{
    /** \brief Sparsities of function inputs and outputs */
    Sparsity get_sparsity_in(casadi_int i) override;
    Sparsity get_sparsity_out(casadi_int i) override;
    /// @}

    ///@{
    /** \brief Names of function input and outputs */
    std::string get_name_in(casadi_int i) override { return simulator_in(i);}
    std::string get_name_out(casadi_int i) override { return simulator_out(i);}
    /// @}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief  Initialize */
    void init(const Dict& opts) override;

    /** Helper for a more powerful 'simulator' factory */
    virtual Function create_advanced(const Dict& opts);

    /** \brief Reset the forward problem */
    virtual void reset(SimulatorMemory* mem, double t,
                       const double* x, const double* z, const double* p) const = 0;

    /** \brief  Advance solution in time */
    virtual void advance(SimulatorMemory* mem, double t,
                         double* x, double* z, double* q) const = 0;

    /** \brief Reset the backward problem */
    virtual void resetB(SimulatorMemory* mem, double t,
                        const double* rx, const double* rz, const double* rp) const = 0;

    /** \brief  Retreat solution in time */
    virtual void retreat(SimulatorMemory* mem, double t,
                         double* rx, double* rz, double* rq) const = 0;

    /** \brief  evaluate */
    int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

    /** \brief  Print solver statistics */
    virtual void print_stats(SimulatorMemory* mem) const {}

    /** \brief  Set stop time for the integration */
    virtual void setStopTime(SimulatorMemory* mem, double tf) const;

    /// Create sparsity pattern of the extended Jacobian (forward problem)
    Sparsity sp_jac_dae();

    /// Create sparsity pattern of the extended Jacobian (backward problem)
    Sparsity sp_jac_rdae();

    // Sparsity pattern of the extended Jacobians
    Sparsity sp_jac_dae_, sp_jac_rdae_;

    ///@{
    // Shorthands
    const Sparsity&  t() const { return oracle_.sparsity_in(DE2_T);}
    const Sparsity&  x() const { return oracle_.sparsity_in(DE2_X);}
    const Sparsity&  z() const { return oracle_.sparsity_in(DE2_Z);}
    const Sparsity&  p() const { return oracle_.sparsity_in(DE2_P);}
    const Sparsity&  q() const { return oracle_.sparsity_out(DE2_QUAD);}
    const Sparsity& rx() const { return oracle_.sparsity_in(DE2_RX);}
    const Sparsity& rz() const { return oracle_.sparsity_in(DE2_RZ);}
    const Sparsity& rp() const { return oracle_.sparsity_in(DE2_RP);}
    const Sparsity& rq() const { return oracle_.sparsity_out(DE2_RQUAD);}
    ///@}

    /// Number of states for the forward integration
    casadi_int nx_, nz_, nq_, nx1_, nz1_, nq1_;

    /// Number of states for the backward integration
    casadi_int nrx_, nrz_, nrq_, nrx1_, nrz1_, nrq1_;

    /// Number of forward and backward parameters
    casadi_int np_, nrp_, np1_, nrp1_;

    /// Number of sensitivities
    casadi_int ns_;

    // Time grid
    std::vector<double> grid_;
    casadi_int ngrid_;

    // Copy of the options
    Dict opts_;

    /// Options
    bool print_stats_;

    /// Output the state at the initial time
    bool output_t0_;
    casadi_int ntout_;

    // Creator function for internal class
    typedef Simulator* (*Creator)(const std::string& name, const Function& oracle);

    // No static functions exposed
    struct Exposed{ };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Infix
    static const std::string infix_;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_SIMULATOR_IMPL_HPP
