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

#ifndef CASADI_ORACLE_HPP
#define CASADI_ORACLE_HPP

#include "function.hpp"

namespace casadi {

#ifndef SWIG
  // A linear combination of inputs
  typedef std::pair<std::string, std::vector<std::string> > LinComb;

  /** Symbolic representation of a problem */
  template<typename XType>
  class CASADI_EXPORT Problem {
  public:
    std::vector<XType> in;
    std::vector<XType> out;
    std::vector<std::string> ischeme;
    std::vector<std::string> oscheme;

    // Input sparsity
    const Sparsity& sparsity_in(int i) { return in.at(i).sparsity();}

    // Output sparsity
    const Sparsity& sparsity_out(int i) { return out.at(i).sparsity();}

    // Constructor (expressions given)
    Problem(const std::vector<XType>& in, const std::vector<XType>& out,
            const std::vector<std::string>& ischeme,
            const std::vector<std::string>& oscheme)
      : in(in), out(out), ischeme(ischeme), oscheme(oscheme) {
    }

    // Constructor (only schemes)
    Problem(const std::vector<std::string>& ischeme,
            const std::vector<std::string>& oscheme)
      : in(ischeme.size()), out(oscheme.size()), ischeme(ischeme), oscheme(oscheme) {
    }

    // Factory
    Function create(const std::string& fname,
                    const std::vector<std::string>& s_in,
                    const std::vector<std::string>& s_out,
                    const std::vector<LinComb>& lincomb,
                    const Dict& opts) const;

    /** \brief Which variables enter nonlinearly */
    std::vector<bool> nl_var(const std::string& s_in,
                             const std::vector<std::string>& s_out);
  };

  typedef Problem<SX> SXProblem;
  typedef Problem<MX> MXProblem;

  /** Can be either an SXProblem or MXProblem */
  class CASADI_EXPORT XProblem {
  public:
    union {
      SXProblem *sx_p;
      MXProblem *mx_p;
    };
    bool is_sx;
    /// Object is SX
    XProblem(const SXProblem& d);
    /// Object is MX
    XProblem(const MXProblem& d);
    /// Copy constructor
    XProblem(const XProblem& d);
    XProblem& operator=(const XProblem& d);
    /// Destructor
    ~XProblem();
    // Input sparsity
    const Sparsity& sparsity_in(int i) const;

    // Output sparsity
    const Sparsity& sparsity_out(int i) const;

    // Factory function
    Function create(const std::string& fname,
                    const std::vector<std::string>& s_in,
                    const std::vector<std::string>& s_out,
                    const std::vector<LinComb>& lincomb = std::vector<LinComb>(),
                    const Dict& opts=Dict()) const;

    /** \brief Which variables enter nonlinearly */
    std::vector<bool> nl_var(const std::string& s_in,
                             const std::vector<std::string>& s_out);
  };
#endif // SWIG

} // namespace casadi

#endif // CASADI_ORACLE_HPP
