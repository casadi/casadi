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


#ifndef CASADI_NLP_BUILDER_HPP
#define CASADI_NLP_BUILDER_HPP

#include "../mx/mx.hpp"

namespace casadi {

  /** \brief A symbolic NLP representation
  \date 2012-2015
  \author Joel Andersson
  */
  class CASADI_EXPORT NlpBuilder : public PrintableObject<NlpBuilder> {
  public:

    /** @name Symbolic representation of the NLP
    *  Data members
    */
    ///@{

    /// Variables
    std::vector<MX> x;

    /// Objective
    MX f;

    /// Constraints
    std::vector<MX> g;

    /// Bounds on x
    std::vector<double> x_lb, x_ub;

    /// Bounds on g
    std::vector<double> g_lb, g_ub;

    /// Primal initial guess
    std::vector<double> x_init;

    /// Dual initial guess
    std::vector<double> lambda_init;

    ///@}

    /// Import an .nl file
    void import_nl(const std::string& filename, const Dict& opts = Dict());

    /// Print a description of the object
    void print(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const;

    /// Print a representation of the object
    void repr(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const;
  };

#ifndef SWIG
  /** \Helper class for .nl import
  The .nl format is described in "Writing .nl Files" paper by David M. Gay (2005)
  \date 2016
  \author Joel Andersson
  */
  class CASADI_EXPORT NlImporter {
  public:
    // Constructor
    NlImporter(NlpBuilder& nlp, const std::string& filename, const Dict& opts);
    // Destructor
    ~NlImporter();
  private:
    // Reference to the class
    NlpBuilder& nlp_;
    // Options
    bool verbose_;
    // File stream
    std::ifstream s_;
    // All variables, including dependent
    std::vector<MX> v_;
    // Number of objectives and constraints
    int n_var_, n_con_, n_obj_, n_eq_, n_lcon_;
    // Parse the file
    void parse();
    // Imported function description
    void F_segment();
    // Suffix values
    void S_segment();
    // Defined variable definition
    void V_segment();
    // Algebraic constraint body
    void C_segment();
    // Logical constraint expression
    void L_segment();
    // Objective function
    void O_segment();
    // Dual initial guess
    void d_segment();
    // Primal initial guess
    void x_segment();
    // Bounds on algebraic constraint bodies ("ranges")
    void r_segment();
    // Bounds on variable
    void b_segment();
    // Jacobian row counts
    void k_segment();
    // Linear terms in the constraint function
    void J_segment();
    // Linear terms in the objective function
    void G_segment();
    /// Read an expression from an NL-file (Polish infix format)
    MX expr();
  };
#endif // SWIG

} // namespace casadi

#endif // CASADI_NLP_BUILDER_HPP
