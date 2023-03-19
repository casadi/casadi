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


#ifndef CASADI_RUNGE_KUTTA_HPP
#define CASADI_RUNGE_KUTTA_HPP

#include "casadi/core/integrator_impl.hpp"
#include <casadi/solvers/casadi_integrator_rk_export.h>

/** \defgroup plugin_Integrator_rk Title
    \par

      Fixed-step explicit Runge-Kutta integrator for ODEs
      Currently implements RK4.

      The method is still under development

    \identifier{23a} */
/** \pluginsection{Integrator,rk} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{Integrator,rk}


  
      @copydoc plugin_Integrator_rk

      \author Joel Andersson
      \date 2011-2014
  */
  class CASADI_INTEGRATOR_RK_EXPORT RungeKutta : public FixedStepIntegrator {
   public:

    /// Constructor
    RungeKutta(const std::string& name, const Function& dae, double t0,
      const std::vector<double>& tout);

    /** \brief  Create a new integrator */
    static Integrator* creator(const std::string& name, const Function& dae,
        double t0, const std::vector<double>& tout) {
      return new RungeKutta(name, dae, t0, tout);
    }

    /// Destructor
    ~RungeKutta() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "rk";}

    // Get name of the class
    std::string class_name() const override { return "RungeKutta";}

    /// Initialize stage
    void init(const Dict& opts) override;

    /// Setup step functions
    void setup_step() override;

    /// A documentation string
    static const std::string meta_doc;

    /// Continuous time dynamics
    Function f_, g_;

    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize into MX */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new RungeKutta(s); }

   protected:

    ///@{
    /** \brief IO conventions for continuous time dynamics */
    enum OdeIn { ODE_T, ODE_X, ODE_P, ODE_U, ODE_NUM_IN};
    enum OdeOut { ODE_ODE, ODE_QUAD, ODE_NUM_OUT};
    enum ROdeIn { RODE_T, RODE_X, RODE_P, RODE_U, RODE_RX, RODE_RP, RODE_NUM_IN};
    enum ROdeOut { RODE_RODE, RODE_RQUAD, RODE_UQUAD, RODE_NUM_OUT};
    ///@}

    /** \brief Deserializing constructor */
    explicit RungeKutta(DeserializingStream& s);
  };

} // namespace casadi

/// \endcond
#endif // CASADI_RUNGE_KUTTA_HPP
