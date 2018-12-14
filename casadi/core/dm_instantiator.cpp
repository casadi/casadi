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

#define CASADI_DM_INSTANTIATOR_CPP
#include "matrix_impl.hpp"

using namespace std;

namespace casadi {


  template<>
  DM CASADI_EXPORT DM::
  solve(const DM& A, const DM& b,
        const string& lsolver, const Dict& dict) {
    Linsol mysolver("tmp", lsolver, A.sparsity(), dict);
    return mysolver.solve(A, b, false);
  }

  template<>
  DM CASADI_EXPORT DM::
  inv(const DM& A,
        const string& lsolver, const Dict& dict) {
    return solve(A, DM::eye(A.size1()), lsolver, dict);
  }

  template<>
  DM CASADI_EXPORT DM::
  pinv(const DM& A, const string& lsolver,
       const Dict& dict) {
    if (A.size1()>=A.size2()) {
      return solve(mtimes(A.T(), A), A.T(), lsolver, dict);
    } else {
      return solve(mtimes(A, A.T()), A, lsolver, dict).T();
    }
  }

  template<>
  DM CASADI_EXPORT DM::
  rand(const Sparsity& sp) { // NOLINT(runtime/threadsafe_fn)
    // C++11 random number generator
    std::uniform_real_distribution<double> distribution(0., 1.);
    // Nonzeros
    std::vector<double> nz(sp.nnz());
    for (double& e : nz) e = distribution(rng_);
    // Construct return object
    return DM(sp, nz, false);
  }

  template<>
  DM CASADI_EXPORT DM::
  expm(const DM& A) {
    Function ret = expmsol("mysolver", "slicot", A.sparsity());
    return ret(std::vector<DM>{A, 1})[0];
  }

  template<>
  DM CASADI_EXPORT DM::
  expm_const(const DM& A, const DM& t) {
    return expm(A*t);
  }

  template<> void CASADI_EXPORT DM::export_code(const std::string& lang,
       std::ostream &stream, const Dict& options) const {

    casadi_assert(lang=="matlab", "Only matlab language supported for now.");

    // Default values for options
    bool opt_inline = false;
    std::string name = "m";
    casadi_int indent_level = 0;
    bool spoof_zero = false;

    // Read options
    for (auto&& op : options) {
      if (op.first=="inline") {
        opt_inline = op.second;
      } else if (op.first=="name") {
        name = op.second.to_string();
      } else if (op.first=="indent_level") {
        indent_level = op.second;
      } else if (op.first=="spoof_zero") {
        spoof_zero = op.second;
      } else {
        casadi_error("Unknown option '" + op.first + "'.");
      }
    }

    // Construct indent string
    std::string indent;
    for (casadi_int i=0;i<indent_level;++i) {
      indent += "  ";
    }

    casadi_assert(!opt_inline, "Inline not supported for now.");

    // Prepare stream for emitting full precision
    std::ios_base::fmtflags fmtfl = stream.flags();
    stream << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1);

    // Obtain nonzeros of matrix
    std::vector<double> d = nonzeros();

    // Spoof numericals
    if (spoof_zero) {
      for (double& e : d) {
        if (e==0) e=1e-200;
      }
    }

    // Short-circuit for (dense) scalars
    if (is_scalar(true)) {
      stream << indent << name << " = " << d[0] << ";" << std::endl;
      stream.flags(fmtfl);
      return;
    }

    // Are all nonzeros equal?
    bool all_equal = true;
    for (double e : d) {
      if (e!=d[0]) {
        all_equal = false;
        break;
      }
    }

    if (all_equal && !d.empty()) {
      // No need to export all individual nonzeros if they are all equal
      stream << indent << name << "_nz = ones(1, " << d.size() << ")*" << d[0] << ";" << std::endl;
    } else {
      // Export nonzeros
      stream << indent << name << "_nz = [";
      for (casadi_int i=0;i<d.size();++i) {
        stream << d[i] << " ";
        if ((i+1)%20 == 0) stream << "..." << std::endl << indent << "  ";
      }
      stream << "];" << std::endl;
    }

    // Reset stream properties
    stream.flags(fmtfl);

    // Cast nonzeros in correct shape
    if (is_dense()) {
      // Special case for dense (for readibility of exported code)
      stream << indent << name << " = reshape(";
      stream << name << "_nz, ";
      stream << size1() << ", " << size2() << ");" << endl;
    } else {
      // For sparse matrices, export Sparsity and use sparse constructor
      Dict opts;
      opts["as_matrix"] = false;
      opts["indent_level"] = indent_level;
      opts["name"] = name;
      opts["indent_level"] = opt_inline;
      sparsity().export_code(lang, stream, opts);
      stream << indent << name << " = sparse(" << name << "_i, " << name << "_j, ";
      stream << name << "_nz, ";
      stream << size1() << ", " << size2() << ");" << endl;
    }
  }

  template<>
  Dict CASADI_EXPORT DM::info() const {
    return {{"sparsity", sparsity().info()}, {"data", nonzeros()}};
  }

  template<>
  void CASADI_EXPORT DM::to_file(const std::string& filename,
      const std::string& format_hint) const {
    std::string format = Sparsity::file_format(filename, format_hint);
    std::ofstream out(filename);
    if (format=="mtx") {
      out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1);
      out << "%%MatrixMarket matrix coordinate real general" << std::endl;
      out << size1() << " " << size2() << " " << nnz() << std::endl;
      std::vector<casadi_int> row = sparsity().get_row();
      std::vector<casadi_int> col = sparsity().get_col();

      for (casadi_int k=0;k<row.size();++k) {
        out << row[k]+1 << " " << col[k]+1 << " " << nonzeros_[k] << std::endl;
      }
    } else {
      casadi_error("Unknown format '" + format + "'");
    }
  }

  // Instantiate templates
  template class CASADI_EXPORT casadi_limits<double>;
  template class CASADI_EXPORT Matrix<double>;

} // namespace casadi
