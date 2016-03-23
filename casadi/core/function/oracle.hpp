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
#include "compiler.hpp"

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

    // Construct from a DLL or C code
    static Oracle* construct(const std::string& fname,
                             const std::string& all_io="all_io");

    // Construct using a compiler instance
    static Oracle* construct(const Compiler& compiler,
                             const std::string& all_io="all_io");

    // Destructor
    virtual ~Oracle() {}

    // Clone
    virtual Oracle* clone() const = 0;

    // Expand MX -> SX
    virtual Oracle* expand() const;

    // Input sparsity
    virtual const Sparsity& sparsity_in(int i) const = 0;

    // Output sparsity
    virtual const Sparsity& sparsity_out(int i) const = 0;

    // Input names
    virtual std::string name_in(int i) const;

    // Output names
    virtual std::string name_out(int i) const;

    // Factory
    virtual Function create(const std::string& fname,
                            const std::vector<std::string>& s_in,
                            const std::vector<std::string>& s_out,
                            const std::vector<LinComb>& lincomb=std::vector<LinComb>(),
                            const Dict& opts=Dict()) const = 0;

    /** \brief Get type name */
    virtual std::string type_name() const = 0;

    /** \brief Which variables enter nonlinearly */
    virtual std::vector<bool> nl_var(const std::string& s_in,
                                     const std::vector<std::string>& s_out) const;

    /** \brief Generate a function with all inputs and outputs */
    virtual Function all_io(const std::string& fname="all_io",
                            const Dict& opts=Dict()) const;
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

    // Input names
    virtual std::string name_in(int i) const { return ischeme_.at(i);}

    // Output names
    virtual std::string name_out(int i) const { return oscheme_.at(i);}

    // Factory
    virtual Function create(const std::string& fname,
                            const std::vector<std::string>& s_in,
                            const std::vector<std::string>& s_out,
                            const std::vector<LinComb>& lincomb,
                            const Dict& opts) const;

    /** \brief Get type name */
    virtual std::string type_name() const;

    /** \brief Which variables enter nonlinearly */
    virtual std::vector<bool> nl_var(const std::string& s_in,
                                     const std::vector<std::string>& s_out) const;

    /** \brief Generate a function with all inputs and outputs */
    virtual Function all_io(const std::string& fname, const Dict& opts) const;
  };

  /** Problem stored externally */
  template<typename LibType>
  class CASADI_EXPORT LibOracle : public Oracle {
  private:
    LibType libtype_;
    Function all_io_;
  public:
    // Constructor
    LibOracle(const LibType& libtype, const std::string& all_io);

    // Destructor
    virtual ~LibOracle() {}

    // Clone
    virtual Oracle* clone() const { return new LibOracle(*this);}

    // Input sparsity
    virtual const Sparsity& sparsity_in(int i) const;

    // Output sparsity
    virtual const Sparsity& sparsity_out(int i) const;

    // Factory
    virtual Function create(const std::string& fname,
                            const std::vector<std::string>& s_in,
                            const std::vector<std::string>& s_out,
                            const std::vector<LinComb>& lincomb,
                            const Dict& opts) const;

    /** \brief Get type name */
    virtual std::string type_name() const;
  };

#endif // SWIG

} // namespace casadi

#endif // CASADI_ORACLE_HPP
