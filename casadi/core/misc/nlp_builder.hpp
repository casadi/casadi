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

#include "../sx/sx_elem.hpp"

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
      SX x;

      /// Objective functions
      SX f;

      /// Constraint functions
      SX g;

      /// Bounds on x
      DM x_lb, x_ub;

      /// Bounds on g
      DM g_lb, g_ub;

      /// Primal initial guess
      DM x_init;

      /// Dual initial guess
      DM lambda_init;

    ///@}

    /// Parse an AMPL och PyOmo NL-file
    void parse_nl(const std::string& filename, const Dict& options = Dict());

    /// Print a description of the object
    void print(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const;

    /// Print a representation of the object
    void repr(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const;

#ifndef SWIG
  protected:

    /// Read an expression from an NL-file (Polish infix format)
    static SXElem read_expr(std::istream &stream, const std::vector<SXElem>& v);

#endif // SWIG
};

} // namespace casadi

#endif // CASADI_NLP_BUILDER_HPP
