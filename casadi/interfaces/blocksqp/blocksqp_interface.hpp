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


#ifndef CASADI_BLOCKSQP_INTERFACE_HPP
#define CASADI_BLOCKSQP_INTERFACE_HPP

#include <casadi/interfaces/blocksqp/casadi_nlpsol_blocksqp_export.h>
#include <blocksqp.hpp>
#include "casadi/core/function/nlpsol_impl.hpp"

/** \defgroup plugin_Nlpsol_blocksqp
  BLOCKSQP interface
*/

/** \pluginsection{Nlpsol,blocksqp} */

/// \cond INTERNAL
namespace casadi {
  // Forward declaration
  class BlocksqpInterface;

  struct CASADI_NLPSOL_BLOCKSQP_EXPORT BlocksqpMemory : public NlpsolMemory {
  };

  // Problem class
  class BlocksqpProblem : public blocksqp::Problemspec {
  public:
    const BlocksqpInterface& self;
    BlocksqpMemory* m;
    BlocksqpProblem(const BlocksqpInterface& self, BlocksqpMemory* m);

    // Set initial values for xi (and possibly lambda) and parts of the
    // Jacobian that correspond to linear constraints (dense version).
    virtual void initialize(blocksqp::Matrix &xi,
                            blocksqp::Matrix &lambda,
                            blocksqp::Matrix &constrJac);

    // Set initial values for xi (and possibly lambda) and parts of the
    // Jacobian that correspond to linear constraints (sparse version).
    virtual void initialize(blocksqp::Matrix &xi,
                            blocksqp::Matrix &lambda,
                            double *&jacNz,
                            int *&jacIndRow,
                            int *&jacIndCol);

    /// Evaluate objective, constraints, and derivatives (dense version).
    virtual void evaluate(const blocksqp::Matrix &xi,
                          const blocksqp::Matrix &lambda,
                          double *objval,
                          blocksqp::Matrix &constr,
                          blocksqp::Matrix &gradObj,
                          blocksqp::Matrix &constrJac,
                          blocksqp::SymMatrix *&hess,
                          int dmode,
                          int *info);

    /// Evaluate objective, constraints, and derivatives (sparse version).
    virtual void evaluate(const blocksqp::Matrix &xi,
                          const blocksqp::Matrix &lambda,
                          double *objval,
                          blocksqp::Matrix &constr,
                          blocksqp::Matrix &gradObj,
                          double *&jacNz,
                          int *&jacIndRow,
                          int *&jacIndCol,
                          blocksqp::SymMatrix *&hess,
                          int dmode,
                          int *info);

    // Generic method to convert dense constraint Jacobian to a sparse matrix
    // in Harwell--Boeing (column compressed) format.
    virtual void convertJacobian(const blocksqp::Matrix &constrJac,
                                 double *&jacNz,
                                 int *&jacIndRow,
                                 int *&jacIndCol,
                                 bool firstCall = 0);
  };


  /** \brief \pluginbrief{Nlpsol,blocksqp}
     @copydoc Nlpsol_doc
     @copydoc plugin_Nlpsol_blocksqp
  */
  class CASADI_NLPSOL_BLOCKSQP_EXPORT BlocksqpInterface : public Nlpsol {
  public:
    explicit BlocksqpInterface(const std::string& name, const Function& nlp);
    virtual ~BlocksqpInterface();

    // Get name of the plugin
    virtual const char* plugin_name() const { return "blocksqp";}

    /** \brief  Create a new NLP Solver */
    static Nlpsol* creator(const std::string& name, const Function& nlp) {
      return new BlocksqpInterface(name, nlp);
    }

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    // Initialize the solver
    virtual void init(const Dict& opts);

    /** \brief Create memory block */
    virtual void* alloc_memory() const { return new BlocksqpMemory();}

    /** \brief Free memory block */
    virtual void free_memory(void *mem) const { delete static_cast<BlocksqpMemory*>(mem);}

    /** \brief Initalize memory block */
    virtual void init_memory(void* mem) const;

    /** \brief Set the (persistent) work vectors */
    virtual void set_work(void* mem, const double**& arg, double**& res,
                          int*& iw, double*& w) const;

    // Solve the NLP
    virtual void solve(void* mem) const;

    /// A documentation string
    static const std::string meta_doc;

    // Block partitioning
    std::vector<int> blocks_;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_BLOCKSQP_INTERFACE_HPP
