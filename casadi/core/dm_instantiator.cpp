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
      const Sparsity& sp, const double* nonzeros,
      const std::string& format_hint) {
    std::string format = Sparsity::file_format(filename, format_hint, {"mtx", "txt"});
    std::ofstream out(filename);
    if (format=="mtx") {
      normalized_setup(out);
      out << "%%MatrixMarket matrix coordinate real general" << std::endl;
      out << sp.size1() << " " << sp.size2() << " " << sp.nnz() << std::endl;
      std::vector<casadi_int> row = sp.get_row();
      std::vector<casadi_int> col = sp.get_col();

      for (casadi_int k=0;k<row.size();++k) {
        out << row[k]+1 << " " << col[k]+1 << " ";
        normalized_out(out, nonzeros ? nonzeros[k]: 0);
        out << std::endl;
      }
    } else if (format=="txt") {
      normalized_setup(out);
      out << std::left;
      // Access data structures
      casadi_int size1 = sp.size1();
      casadi_int size2 = sp.size2();
      const casadi_int* colind = sp.colind();
      const casadi_int* row = sp.row();

      // Index counter for each column
      std::vector<casadi_int> ind(colind, colind+size2+1);

      // Make enough room to put -3.3e-310
      casadi_int w = std::numeric_limits<double>::digits10 + 9;

      // Loop over rows
      for (casadi_int rr=0; rr<size1; ++rr) {
        // Loop over columns
        for (casadi_int cc=0; cc<size2; ++cc) {
          // Set filler execptfor last column
          if (cc<size2-1) out << std::setw(w);
          // String representation of element
          if (ind[cc]<colind[cc+1] && row[ind[cc]]==rr) {
            normalized_out(out, nonzeros ? nonzeros[ind[cc]++]: 0);
          } else {
            out << std::setw(w) << "00";
          }
          if (cc<size2-1) out << " ";
        }
        out << std::endl;
      }
    } else {
      casadi_error("Unknown format '" + format + "'");
    }
  }

  template<>
  DM CASADI_EXPORT DM::from_file(const std::string& filename, const std::string& format_hint) {
    std::string format = Sparsity::file_format(filename, format_hint, {"mtx", "txt"});
    std::ifstream in(filename);
    casadi_assert(in.good(), "Could not open '" + filename + "'.");

    if (format=="txt") {
      std::string line;
      std::vector<double> values;
      casadi_int n_row = 0;
      casadi_int n_col = 0;
      bool first_line = true;
      std::istringstream stream;

      std::vector<casadi_int> row;
      std::vector<casadi_int> col;

      normalized_setup(stream);

      // Read line-by-line
      while (std::getline(in, line)) {
        // Ignore empty lines
        if (line.empty()) continue;

        // Ignore lines with comments
        if (line[0]=='%' || line[0]=='#' || line[0]=='/') continue;

        // Populate a stream for pulling doubles
        stream.clear();
        stream.str(line);

        // Keep pulling doubles from line
        double val;
        casadi_int i=0;
        for (i=0; !stream.eof(); ++i) {
          casadi_int start = stream.tellg();
          int ret = normalized_in(stream, val);

          if (ret==-1) break; // EOL reached
          casadi_assert(ret==0, "Parsing error on line " + str(i+1) + ", column " + str(start+1));
          casadi_int stop = line.size();
          if (!stream.eof()) stop = stream.tellg();

          // Check if structural zero
          bool structural_zero = false;
          if (val==0) {
            // Check if stream contained '00'
            casadi_int n_zeros = 0;
            for (casadi_int k=start;k<stop;++k) {
              char c = line.at(k);
              if (c==' ' || c=='\t') continue;
              if (c=='0') {
                n_zeros++;
              } else {
                break;
              }
            }
            if (n_zeros==2) structural_zero = true;
          }

          if (!structural_zero) {
            row.push_back(n_row);
            col.push_back(i);
            values.push_back(val);
          }
          if (first_line) n_col++;
        }

        // Dimension check
        casadi_assert(i==n_col, "Inconsistent dimensions. "
        "File started with " + str(n_col) + ", while line " + str(n_row+1) +
        " has " + str(i) + ".");

        first_line = false;
        n_row++;
      }
      return DM::triplet(row, col, values, n_row, n_col);
    } else if (format=="mtx") {
      std::string line;
      bool first_line = true;
      std::istringstream stream;

      casadi_int n_row=0, n_col=0, nnz=0;
      std::vector<double> values;
      std::vector<casadi_int> row, col;
      normalized_setup(stream);

      casadi_int i=0;
      // Read line-by-line
      while (std::getline(in, line)) {
        // Ignore empty lines
        if (line.empty()) continue;

        // Ignore lines with comments
        if (line[0]=='%' || line[0]=='#' || line[0]=='/') continue;

        // Populate a stream for pulling doubles
        stream.clear();
        stream.str(line);

        if (first_line) {
          stream >> n_row;
          stream >> n_col;
          stream >> nnz;
          casadi_assert(!stream.fail(), "Could not parse first line");
          values.reserve(nnz);
          row.reserve(nnz);
          col.reserve(nnz);
          first_line = false;

        } else {
          casadi_int r, c;
          double val;
          stream >> r;
          stream >> c;
          casadi_assert(normalized_in(stream, val)==0, "Parse error");
          row.push_back(r-1);
          col.push_back(c-1);
          values.push_back(val);
          i++;
        }
      }
      return DM::triplet(row, col, values, n_row, n_col);
    } else {
      casadi_error("Unknown format '" + format + "'");
    }
  }


  // Instantiate templates
  template class CASADI_EXPORT casadi_limits<double>;
  template class CASADI_EXPORT Matrix<double>;

} // namespace casadi
