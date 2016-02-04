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

  /** Oracle abstract base class */
  class CASADI_EXPORT Oracle {
    public:
    // Construct from SX expression
    static Oracle* construct(const std::vector<SX>& in,
                             const std::vector<SX>& out,
                             const std::vector<std::string>& ischeme,
                             const std::vector<std::string>& oscheme);

    // Construct from MX expression
    static Oracle* construct(const std::vector<MX>& in,
                             const std::vector<MX>& out,
                             const std::vector<std::string>& ischeme,
                             const std::vector<std::string>& oscheme);

    // Destructor
    virtual ~Oracle() {}

    // Clone
    virtual Oracle* clone() const = 0;

    // Input sparsity
    virtual const Sparsity& sparsity_in(int i) const = 0;

    // Output sparsity
    virtual const Sparsity& sparsity_out(int i) const = 0;

    // Factory
    virtual Function create(const std::string& fname,
                            const std::vector<std::string>& s_in,
                            const std::vector<std::string>& s_out,
                            const std::vector<LinComb>& lincomb={},
                            const Dict& opts=Dict()) const = 0;

    /** \brief Which variables enter nonlinearly */
    virtual std::vector<bool> nl_var(const std::string& s_in,
                                     const std::vector<std::string>& s_out) = 0;
  };

  /** Symbolic representation of a problem */
  template<typename XType>
  class CASADI_EXPORT XOracle : public Oracle {
  private:
    std::vector<XType> in_;
    std::vector<XType> out_;
    std::vector<std::string> ischeme_;
    std::vector<std::string> oscheme_;
  public:
    // Constructor
    XOracle(const std::vector<XType>& in, const std::vector<XType>& out,
            const std::vector<std::string>& ischeme,
            const std::vector<std::string>& oscheme)
      : in_(in), out_(out), ischeme_(ischeme), oscheme_(oscheme) {
    }

    // Destructor
    virtual ~XOracle() {}

    // Clone
    virtual Oracle* clone() const { return new XOracle(*this);}

    // Input sparsity
    virtual const Sparsity& sparsity_in(int i) const { return in_.at(i).sparsity();}

    // Output sparsity
    virtual const Sparsity& sparsity_out(int i) const { return out_.at(i).sparsity();}

    // Factory
    virtual Function create(const std::string& fname,
                            const std::vector<std::string>& s_in,
                            const std::vector<std::string>& s_out,
                            const std::vector<LinComb>& lincomb,
                            const Dict& opts) const;

    /** \brief Which variables enter nonlinearly */
    virtual std::vector<bool> nl_var(const std::string& s_in,
                                     const std::vector<std::string>& s_out);
  };

#endif // SWIG

} // namespace casadi

#endif // CASADI_ORACLE_HPP
