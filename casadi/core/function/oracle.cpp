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


#include "oracle.hpp"
#include "external.hpp"

using namespace std;

namespace casadi {

  Oracle* Oracle::construct(const std::vector<SX>& in,
                            const std::vector<SX>& out,
                            const std::vector<std::string>& ischeme,
                            const std::vector<std::string>& oscheme) {
    return new XOracle<SX>(in, out, ischeme, oscheme);
  }

  Oracle* Oracle::construct(const std::vector<MX>& in,
                            const std::vector<MX>& out,
                            const std::vector<std::string>& ischeme,
                            const std::vector<std::string>& oscheme) {
    return new XOracle<MX>(in, out, ischeme, oscheme);
  }

  Oracle* Oracle::construct(const Importer& compiler, const std::string& all_io) {
    return new LibOracle<Importer>(compiler, all_io);
  }

  Oracle* Oracle::construct(const std::string& fname, const std::string& all_io) {
    // If fname ends with .c, JIT
    if (fname.size()>2 && fname.compare(fname.size()-2, fname.size(), ".c")==0) {
      Importer compiler(fname, "clang");
      return construct(compiler, all_io);
    } else {
      return new LibOracle<std::string>(fname, all_io);
    }
  }

  std::string Oracle::name_in(int i) const {
    casadi_error("'name_in' not defined for " + type_name());
  }

  std::string Oracle::name_out(int i) const {
    casadi_error("'name_out' not defined for " + type_name());
  }

  template<typename XType>
  std::string XOracle<XType>::type_name() const {
    return XType::type_name() + "Oracle";
  }

  Function Oracle::all_io(const std::string& fname, const Dict& opts) const {
    casadi_error("'all_io' not defined for " + type_name());
  }

  template<typename XType>
  Function XOracle<XType>::all_io(const std::string& fname, const Dict& opts) const {
    return Function(fname, in_, out_, ischeme_, oscheme_, opts);
  }

  template<typename LibType>
  LibOracle<LibType>::LibOracle(const LibType& libtype, const std::string& all_io)
    : libtype_(libtype) {
    all_io_ = external(all_io, libtype);
  }

  template<typename LibType>
  const Sparsity& LibOracle<LibType>::sparsity_in(int i) const {
    return all_io_.sparsity_in(i);
  }

  template<typename LibType>
  const Sparsity& LibOracle<LibType>::sparsity_out(int i) const {
    return all_io_.sparsity_out(i);
  }

  template<typename LibType>
  std::string LibOracle<LibType>::type_name() const {
    return "LibOracle";
  }

  Oracle* Oracle::expand() const {
    // Get a MX function from all inputs to all outputs
    Function f = all_io();

    // Expand to SX
    if (!f.is_a("sxfunction")) {
      f = f.expand();
    }

    // Get expressions
    vector<SX> arg = f.sx_in();
    vector<SX> res = f(arg);

    // Construct SX oracle
    return construct(arg, res, f.name_in(), f.name_out());
  }

} // namespace casadi
