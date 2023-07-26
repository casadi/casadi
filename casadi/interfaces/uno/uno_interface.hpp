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


#ifndef CASADI_UNO_INTERFACE_HPP
#define CASADI_UNO_INTERFACE_HPP

#include <casadi/interfaces/uno/casadi_nlpsol_uno_export.h>
#include <optimization/Model.hpp>
#include <linear_algebra/RectangularMatrix.hpp>

#include "casadi/core/nlpsol_impl.hpp"

/** \defgroup plugin_Nlpsol_uno Title
    \par

    David Kiessling

  Uno interface

    \identifier{22c} */

/** \pluginsection{Nlpsol,knitro} */



/// \cond INTERNAL
namespace casadi {
  // Forward declaration
  class UnoInterface;
  class CasadiModel;

  /*------------------------
  Definition of UnoMemory
  ------------------------*/

  struct CASADI_NLPSOL_UNO_EXPORT UnoMemory : public NlpsolMemory {
    const UnoInterface& self;

    CasadiModel* model;
    const char* return_status;
    /// Constructor
    UnoMemory(const UnoInterface& uno_interface);

    /// Destructor
    ~UnoMemory();
  };


/* ----------------------------
Definition of Class CasadiModel
-----------------------------*/

class CasadiModel: public Model{

public:
   explicit CasadiModel(const std::string& file_name, const UnoInterface& uno_interface, UnoMemory* mem);
   ~CasadiModel() override {}

   // objective
   [[nodiscard]] double evaluate_objective(const std::vector<double>& x) const override;
   void evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const override;
   // constraints
   void evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const override;
   void evaluate_constraint_gradient(const std::vector<double>& x, size_t j, SparseVector<double>& gradient) const override;
   void evaluate_constraint_jacobian(const std::vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const override;
   // Hessian
   void evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
         SymmetricMatrix<double>& hessian) const override;

   [[nodiscard]] double get_variable_lower_bound(size_t i) const override;
   [[nodiscard]] double get_variable_upper_bound(size_t i) const override;
   [[nodiscard]] BoundType get_variable_bound_type(size_t i) const override;
   [[nodiscard]] double get_constraint_lower_bound(size_t j) const override;
   [[nodiscard]] double get_constraint_upper_bound(size_t j) const override;
   [[nodiscard]] FunctionType get_constraint_type(size_t j) const override;
   [[nodiscard]] BoundType get_constraint_bound_type(size_t j) const override;

   [[nodiscard]] size_t get_number_objective_gradient_nonzeros() const override;
   [[nodiscard]] size_t get_number_jacobian_nonzeros() const override;
   [[nodiscard]] size_t get_number_hessian_nonzeros() const override;

   void get_initial_primal_point(std::vector<double>& x) const override;
   void get_initial_dual_point(std::vector<double>& multipliers) const override;
   void postprocess_solution(Iterate& iterate, TerminationStatus termination_status) const override;

   [[nodiscard]] const std::vector<size_t>& get_linear_constraints() const override;

private:
   UnoMemory* mem_;

   mutable std::vector<double> casadi_tmp_gradient{};
   mutable std::vector<double> casadi_tmp_multipliers{};
   mutable std::vector<double> casadi_tmp_constraint_jacobian{};
   mutable std::vector<double> casadi_tmp_hessian{};

   std::vector<Interval> variables_bounds;
   std::vector<Interval> constraint_bounds;
   std::vector<BoundType> variable_status; /*!< Status of the variables (EQUALITY, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES) */
   std::vector<FunctionType> constraint_type; /*!< Types of the constraints (LINEAR, QUADRATIC, NONLINEAR) */
   std::vector<BoundType> constraint_status; /*!< Status of the constraints (EQUAL_BOUNDS, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES,
 * UNBOUNDED) */

   std::vector<size_t> linear_constraints;

   void generate_variables();
   void generate_constraints();
  //  void set_function_types(std::string file_name);

   void set_number_hessian_nonzeros();
   [[nodiscard]] size_t compute_hessian_number_nonzeros(double objective_multiplier, const std::vector<double>& multipliers) const;
};


  /*------------------------------
  Definition of class UnoInterface
  -------------------------------*/


  /** \brief \pluginbrief{Nlpsol,uno}
     @copydoc Nlpsol_doc
     @copydoc plugin_Nlpsol_uno
  */
  class CASADI_NLPSOL_UNO_EXPORT UnoInterface : public Nlpsol {
  public:

    explicit UnoInterface(const std::string& name, const Function& nlp);
    ~UnoInterface() override;

    // Hessian Sparsity
   Sparsity hesslag_sp_;

   // Jacobian sparsity
   Sparsity jacg_sp_;

    // Get name of the plugin
    const char* plugin_name() const override { return "uno";}

    // Get name of the class
    std::string class_name() const override { return "UnoInterface";}

    /** \brief  Create a new NLP Solver */
    static Nlpsol* creator(const std::string& name, const Function& nlp) {
      return new UnoInterface(name, nlp);
    }

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    // Initialize the solver
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new UnoMemory(*this);}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<UnoMemory*>(mem);}

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    // Solve the NLP
    int solve(void* mem) const override;

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    // UNO options
    Dict opts_;

    /// A documentation string
    static const std::string meta_doc;


  };

} // namespace casadi

/// \endcond
#endif // CASADI_UNO_INTERFACE_HPP
