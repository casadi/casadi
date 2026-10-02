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



#ifndef CASADI_ACADOS_INTEGRATOR_HPP
#define CASADI_ACADOS_INTEGRATOR_HPP

#include <casadi/interfaces/acados/casadi_integrator_acados_export.h>
#include "casadi/core/integrator_impl.hpp"
#include "casadi/core/oracle_function.hpp"
#include "casadi/core/plugin_interface.hpp"

#include <acados_c/sim_interface.h>
#include <acados_c/external_function_interface.h>
#include <acados/sim/sim_erk_integrator.h>
#include <acados/sim/sim_irk_integrator.h>
#include <blasfeo_d_aux.h>

#include <cstring>
#include <memory>

/// \cond INTERNAL
namespace casadi {

  // acados glue, shared with generated code
  #include "acados_sim.hpp"
  #include "acados_chain.hpp"
  #include "acados_runtime.hpp"

  /** \brief Glue between a CasADi DAE and acados sim

      The acados model functions are CasADi Functions registered on an OracleFunction;
      acados calls them through external_function_generic shims (calc_function).
      The acados objects themselves live in work vectors, see acados_sim.hpp.
      Derivatives w.r.t. differentiable p come from separate instances with p in the acados
      controls, w = [u; p], see get_jacobian, get_reverse. Otherwise p never reaches acados.
  */
  class CASADI_INTEGRATOR_ACADOS_EXPORT AcadosModel {
  public:
    void init(const Function& dae, const Dict& opts, bool p_diff);
    // Model function k of variant v, built on first use
    const Function& fcn(bool v, casadi_int k) const;
    // Name of model function k of variant v
    std::string fcn_name(bool v, casadi_int k) const;
    // Register the first n model functions of variant v
    void register_functions(OracleFunction* owner, bool v, casadi_int n) const;
    // sim_collocation_type of the scheme
    int collocation_type_enum() const;
    // Dimensions; parameters in the acados controls of variant 1 (np if p_diff, else 0)
    casadi_int nx, nz, np, nu, npd;
    bool irk, p_diff;
    casadi_int num_stages, num_steps, newton_iter;
    std::string collocation_type;
    Function dae;
    // Model functions: variant 0 with p passed by CasADi, 1 with p in the acados controls
    mutable Function fcns[2][4];
    std::vector<std::string> fields;
    static const std::vector<std::string> option_names;
  };

  /** \brief Integrator plugin that only exists to hand out a AcadosFunction

      integrator(..., 'acados', ...) returns an AcadosFunction from create_advanced,
      like FixedStepIntegrator does for 'simplify'.
  */
  class CASADI_INTEGRATOR_ACADOS_EXPORT AcadosInterface : public Integrator {
  public:
    AcadosInterface(const std::string& name, const Function& dae,
      double t0, const std::vector<double>& tout);

    static Integrator* creator(const std::string& name, const Function& dae,
        double t0, const std::vector<double>& tout) {
      return new AcadosInterface(name, dae, t0, tout);
    }

    ~AcadosInterface() override;

    const char* plugin_name() const override { return "acados";}
    std::string class_name() const override { return "AcadosInterface";}

    static const Options options_;
    const Options& get_options() const override { return options_;}

    Function create_advanced(const Dict& opts) override;

    // Never evaluated: create_advanced returns a AcadosFunction
    int advance_noevent(IntegratorMemory* mem) const override;
    void resetB(IntegratorMemory* mem) const override;
    void impulseB(IntegratorMemory* mem,
      const double* adj_x, const double* adj_z, const double* adj_q) const override;
    void retreat(IntegratorMemory* mem, const double* u,
      double* adj_x, double* adj_p, double* adj_u) const override;

    static const std::string meta_doc;
  };

  struct CASADI_INTEGRATOR_ACADOS_EXPORT AcadosFunctionMemory : public OracleMemory {
    casadi_acados_data<double> d;
  };

  /** \brief acados sim as a Function with the integrator signature

      Derivative strategy of https://github.com/FreyJo/casados-integrators.

      One class, four modes, each a single acados solve:
      NOM: integrator I/O; JAC: jacobian of NOM (S_forw);
      ADJ: reverse of NOM (S_adj, one solve per direction);
      HESS: jacobian of ADJ with nadj=1 (S_hess and S_forw).
      Multiple output times are chained, see acados_chain.hpp.
      Algebraic states are handled by acados (IRK); zf is not supported (nan, acados only
      reports z at the start of an interval).
  */
  class CASADI_INTEGRATOR_ACADOS_EXPORT AcadosFunction
      : public OracleFunction, public PluginInterface<Integrator> {
  public:
    enum Mode {NOM, JAC, ADJ, HESS};

    AcadosFunction(const std::string& name, const Function& dae,
      const std::shared_ptr<AcadosModel>& model, const Dict& model_opts, const Dict& fun_opts,
      const std::vector<double>& T, Mode mode, casadi_int nadj, bool pv,
      const std::vector<Sparsity>& sp_in_nom, const std::vector<Sparsity>& sp_out_nom,
      const std::vector<std::string>& names_in, const std::vector<std::string>& names_out);

    ~AcadosFunction() override;

    const char* plugin_name() const override { return "acados";}
    std::string class_name() const override { return "AcadosFunction";}

    size_t get_n_in() override { return names_in_.size();}
    size_t get_n_out() override { return names_out_.size();}
    std::string get_name_in(casadi_int i) override { return names_in_.at(i);}
    std::string get_name_out(casadi_int i) override { return names_out_.at(i);}
    Sparsity get_sparsity_in(casadi_int i) override;
    Sparsity get_sparsity_out(casadi_int i) override;

    void init(const Dict& opts) override;
    void set_work(void* mem, const double**& arg, double**& res,
      casadi_int*& iw, double*& w) const override;
    void* alloc_mem() const override { return new AcadosFunctionMemory();}
    int init_mem(void* mem) const override;
    void free_mem(void* mem) const override;

    int eval(const double** arg, double** res, casadi_int* iw, double* w,
      void* mem) const override;

    ///@{
    /** \brief Derivatives as in acados: jacobian, reverse, jacobian of reverse */
    bool has_jacobian() const override { return mode_==NOM || (mode_==ADJ && nadj_==1);}
    Function get_jacobian(const std::string& name, const std::vector<std::string>& inames,
      const std::vector<std::string>& onames, const Dict& opts) const override;
    bool has_reverse(casadi_int nadj) const override { return mode_==NOM;}
    Function get_reverse(casadi_int nadj, const std::string& name,
      const std::vector<std::string>& inames, const std::vector<std::string>& onames,
      const Dict& opts) const override;
    ///@}

    ///@{
    /** \brief Code generation */
    bool has_codegen() const override { return true;}
    void codegen_declarations(CodeGenerator& g) const override;
    void codegen_body(CodeGenerator& g) const override;
    ///@}

    ///@{
    /** \brief Serialization, routed through the Integrator plugin registry */
    void serialize_body(SerializingStream &s) const override;
    void serialize_type(SerializingStream &s) const override;
    std::string serialize_base_function() const override { return "Integrator"; }
    static ProtoFunction* deserialize(DeserializingStream& s);
    ///@}

    // Offset of nominal input i in w = [x0; u; p], -1 if no sensitivities
    casadi_int acados_offset(casadi_int i) const;
    // Number of output times
    casadi_int nt() const { return T_.size();}
    // Number of model functions used: the hessian only for HESS
    casadi_int nfun() const { return mode_ == HESS ? 4 : 3;}
    // Fill p_ (except sim_bytes)
    void set_acados_prob();
    // Bytes taken by the acados objects of variant v, measured by laying them out once
    casadi_int measure_sim() const;
    // Derivative function of this family
    Function derivative(const std::string& name, Mode mode, casadi_int nadj, bool pv,
      const std::vector<std::string>& inames, const std::vector<std::string>& onames,
      const Dict& opts) const;
    // Inlined split: blocks w.r.t. p from f1 (p in the acados controls), others from f0
    static Function split_p(const std::string& name, const Function& f0, const Function& f1,
      const std::vector<bool>& from_f1, const std::vector<std::string>& inames,
      const std::vector<std::string>& onames, const Dict& opts);
    // Number of nonzeros of nominal input i
    casadi_int nnz_nom_in(casadi_int i) const { return sp_in_nom_.at(i).nnz();}
    // Index of the adjoint seed on xf among the inputs of ADJ
    static casadi_int adj_seed_xf() {
      return INTEGRATOR_NUM_IN + INTEGRATOR_NUM_OUT + INTEGRATOR_XF;
    }

    std::shared_ptr<AcadosModel> model_;
    Dict model_opts_, fun_opts_;
    // Interval lengths, one per output time
    std::vector<double> T_;
    Mode mode_;
    casadi_int nadj_;
    // p in the acados controls: sensitivities w.r.t. p
    bool pv_;
    std::vector<Sparsity> sp_in_nom_, sp_out_nom_;
    std::vector<std::string> names_in_, names_out_;
    // acados glue problem, model function shapes
    casadi_acados_prob<double> p_;
    std::vector<casadi_int> dims_[4];

  protected:
    explicit AcadosFunction(DeserializingStream& s);
  };

} // namespace casadi
/// \endcond

#endif // CASADI_ACADOS_INTEGRATOR_HPP
