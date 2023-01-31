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


#ifndef CASADI_PROXQP_INTERFACE_HPP
#define CASADI_PROXQP_INTERFACE_HPP

#include "casadi/core/conic_impl.hpp"
#include <casadi/interfaces/proxqp/casadi_conic_proxqp_export.h>

// Proxqp header
#include <proxsuite/proxqp/dense/dense.hpp> // NOLINT(build/include)
#include <proxsuite/proxqp/sparse/sparse.hpp> // NOLINT(build/include)
#include <proxsuite/proxqp/settings.hpp> // NOLINT(build/include)

#include <Eigen/Dense>
#include <Eigen/Sparse>


/** \defgroup plugin_Conic_osqp
    Interface to the PROXQP Solver for quadratic programming
*/

/** \pluginsection{Conic,proxqp} */

/// \cond INTERNAL
namespace casadi {

  struct CASADI_CONIC_PROXQP_EXPORT ProxqpMemory : public ConicMemory {

    typedef Eigen::Triplet<double> T;
    // Solvers
    proxsuite::proxqp::sparse::QP<double, long long> sparse_solver;
    proxsuite::proxqp::dense::QP<double> dense_solver;

    // Working structures:
    Eigen::VectorXd g_vector;    // linear term
    Eigen::VectorXd b_vector;    // equality constraints

    Eigen::VectorXd uba_vector;  // upper bound inequality constraints for Ax
    Eigen::VectorXd lba_vector;  // lower bound inequality constraints for Ax
    Eigen::VectorXd ubx_vector;  // upper bound inequality constraints on variable x
    Eigen::VectorXd lbx_vector;  // upper bound inequality constraints on variable x
    Eigen::VectorXd ub_vector;   // upper bound stacked (Ax, x) input for proxqp
    Eigen::VectorXd lb_vector;   // lower bound stacked (Ax, x) input for proxqp

    // For conversion from casadi::Sparsity to Eigen::SparseMatrix
    std::vector<T> tripletList;  // Used for H and C
    std::vector<T> tripletListEq;  // Used for A
    std::vector<casadi_int> row, col;

    // Results
    std::unique_ptr<Eigen::VectorXd> results_x;
    std::unique_ptr<Eigen::VectorXd> results_y;
    std::unique_ptr<Eigen::VectorXd> results_z;
    double objValue;
    proxsuite::proxqp::QPSolverOutput status;

    /// Constructor
    ProxqpMemory();

    /// Destructor
    ~ProxqpMemory();
  };

  /** \brief \pluginbrief{Conic,proxqp}

      @copydoc Conic_doc
      @copydoc plugin_Conic_osqp

  */
  class CASADI_CONIC_PROXQP_EXPORT ProxqpInterface : public Conic {
  public:
    /** \brief  Create a new Solver */
    explicit ProxqpInterface(const std::string& name,
                             const std::map<std::string, Sparsity>& st);

    /** \brief  Create a new QP Solver */
    static Conic* creator(const std::string& name,
                                     const std::map<std::string, Sparsity>& st) {
      return new ProxqpInterface(name, st);
    }

    /** \brief  Destructor */
    ~ProxqpInterface() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "proxqp";}

    // Get name of the class
    std::string class_name() const override { return "ProxqpInterface";}

    ///@{
    /** \brief const Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief  Initialize */
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new ProxqpMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<ProxqpMemory*>(mem);}

    /// Solve the QP
    int solve(const double** arg, double** res,
        casadi_int* iw, double* w, void* mem) const override;

    /// Can discrete variables be treated
    bool integer_support() const override { return false;}

    /// Can psd constraints be treated
    bool psd_support() const override { return false;}

    /// A documentation string
    static const std::string meta_doc;

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    proxsuite::proxqp::Settings<double> settings_;

    bool warm_start_primal_, warm_start_dual_;
    bool sparse_backend;
    double max_iter;

    // Number of nonzeros in Hessian
    casadi_int nH_;

    // Number of nonzeros in constraint matrix
    casadi_int nA_;

    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize with type disambiguation */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new ProxqpInterface(s); }

  protected:
     /** \brief Deserializing constructor */
    explicit ProxqpInterface(DeserializingStream& e);
  };

} // namespace casadi

/// \endcond
#endif // CASADI_PROXQP_INTERFACE_HPP
