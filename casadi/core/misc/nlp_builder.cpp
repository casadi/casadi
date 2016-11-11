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


#include "nlp_builder.hpp"
#include "../core.hpp"
#include <fstream>

using namespace std;
namespace casadi {

  void NlpBuilder::import_nl(const std::string& filename, const Dict& opts) {
    // Redirect to helper class
    NlImporter(*this, filename, opts);
  }

  void NlpBuilder::print(std::ostream &stream, bool trailing_newline) const {
    stream << "NLP:" << endl;
    stream << "x = " << this->x << endl;
    stream << "f = " << this->f << endl;
    stream << "g = " << this->g << endl;
    if (trailing_newline) stream << endl;
  }

  void NlpBuilder::repr(std::ostream &stream, bool trailing_newline) const {
    stream << "NLP(#x=" << this->x.size() << ", #g=" << this->g.size() << ")";
    if (trailing_newline) stream << endl;
  }

  NlImporter::NlImporter(NlpBuilder& nlp, const std::string& filename, const Dict& opts)
  : nlp_(nlp) {
    // Set default options
    verbose_=false;

    // Read user options
    for (auto&& op : opts) {
      if (op.first == "verbose") {
        verbose_ = op.second;
      } else {
        stringstream ss;
        ss << "Unknown option \"" << op.first << "\"" << endl;
        throw CasadiException(ss.str());
      }
    }
    // Open file for reading
    s_.open(filename.c_str());
    if (verbose_) userOut() << "Reading file \"" << filename << "\"" << endl;

    // Read the header of the NL-file (first 10 lines)
    const int header_sz = 10;
    vector<string> header(header_sz);
    for (int k=0; k<header_sz; ++k) {
      getline(s_, header[k]);
    }

    // Assert that the file is not in binary form
    casadi_assert_message(header.at(0).at(0)=='g',
    "File could not be read, or file is binary format "
    "(currently not supported)");

    // Get the number of objectives and constraints
    stringstream ss(header[1]);
    ss >> n_var_ >> n_con_ >> n_obj_ >> n_eq_ >> n_lcon_;
    if (verbose_) {
      userOut() << "n_var = " << n_var_ << ", n_con  = " << n_con_ << ", n_obj = " << n_obj_
      << ", n_eq = " << n_eq_ << ", n_lcon = " << n_lcon_ << endl;
    }

    // Allocate variables
    nlp_.x = MX::sym("x", 1, 1, n_var_);

    // Allocate f and c
    nlp_.f = 0;
    nlp_.g.resize(n_con_, 0);

    // Allocate bounds for x and primal initial guess
    nlp_.x_lb.resize(n_var_, -inf);
    nlp_.x_ub.resize(n_var_,  inf);
    nlp_.x_init.resize(n_var_, 0);

    // Allocate bounds for g and dual initial guess
    nlp_.g_lb.resize(n_con_, -inf);
    nlp_.g_ub.resize(n_con_,  inf);
    nlp_.lambda_init.resize(n_con_, 0);

    // All variables, including dependent
    v_ = nlp_.x;

    // Read segments
    parse();
  }

  NlImporter::~NlImporter() {
    // Close the NL file
    s_.close();
  }

  void NlImporter::parse() {
    // Segment key
    char key;

    // Process segments
    while (true) {
      // Read segment key
      s_ >> key;
      if (s_.eof()) break; // end of file encountered
      switch (key) {
        case 'F': F_segment(); break;
        case 'S': S_segment(); break;
        case 'V': V_segment(); break;
        case 'C': C_segment(); break;
        case 'L': L_segment(); break;
        case 'O': O_segment(); break;
        case 'd': d_segment(); break;
        case 'x': x_segment(); break;
        case 'r': r_segment(); break;
        case 'b': b_segment(); break;
        case 'k': k_segment(); break;
        case 'J': J_segment(); break;
        case 'G': G_segment(); break;
        default: casadi_error("Unknown .nl segment");
      }
    }
  }

  MX NlImporter::expr() {
    // Read the instruction
    char inst;
    s_ >> inst;

    // Temporaries
    int i;
    double d;

    // Error message
    stringstream msg;

    // Process instruction
    switch (inst) {
      // Symbolic variable
      case 'v':
      // Read the variable number
      s_ >> i;

      // Return the corresponding expression
      return v_.at(i);

      // Numeric expression
      case 'n':

      // Read the floating point number
      s_ >> d;

      // Return an expression containing the number
      return d;

      // Operation
      case 'o':

      // Read the operation
      s_ >> i;

      // Process
      switch (i) {

        // Unary operations, class 1 in Gay2005
        case 13:  case 14:  case 15:  case 16:  case 34:  case 37:  case 38:  case 39:  case 40:
        case 41:  case 43:  case 42:  case 44:  case 45:  case 46:  case 47:  case 49:  case 50:
        case 51:  case 52:  case 53:
        {
          // Read dependency
          MX x = expr();

          // Perform operation
          switch (i) {
            case 13:  return floor(x);
            case 14:  return ceil(x);
            case 15:  return abs(x);
            case 16:  return -x;
            case 34:  return logic_not(x);
            case 37:  return tanh(x);
            case 38:  return tan(x);
            case 39:  return sqrt(x);
            case 40:  return sinh(x);
            case 41:  return sin(x);
            case 42:  return log10(x);
            case 43:  return log(x);
            case 44:  return exp(x);
            case 45:  return cosh(x);
            case 46:  return cos(x);
            // case 47:  return atanh(x); FIXME
            case 49:  return atan(x);
            // case 50:  return asinh(x); FIXME
            case 51:  return asin(x);
            // case 52:  return acosh(x); FIXME
            case 53:  return acos(x);

            default:
            msg << "Unknown unary operation: \"" << i << "\"";
          }
          break;
        }

        // Binary operations, class 2 in Gay2005
        case 0:   case 1:   case 2:   case 3:   case 4:   case 5:   case 6:   case 20:  case 21:
        case 22:  case 23:  case 24:  case 28:  case 29:  case 30:  case 48:  case 55:  case 56:
        case 57:  case 58:  case 73:
        {
          // Read dependencies
          MX x = expr();
          MX y = expr();

          // Perform operation
          switch (i) {
            case 0:   return x + y;
            case 1:   return x - y;
            case 2:   return x * y;
            case 3:   return x / y;
            // case 4:   return rem(x, y); FIXME
            case 5:   return pow(x, y);
            // case 6:   return x < y; // TODO(Joel): Verify this,
            // what is the difference to 'le' == 23 below?
            case 20:  return logic_or(x, y);
            case 21:  return logic_and(x, y);
            case 22:  return x < y;
            case 23:  return x <= y;
            case 24:  return x == y;
            case 28:  return x >= y;
            case 29:  return x > y;
            case 30:  return x != y;
            case 48:  return atan2(x, y);
            // case 55:  return intdiv(x, y); // FIXME
            // case 56:  return precision(x, y); // FIXME
            // case 57:  return round(x, y); // FIXME
            // case 58:  return trunc(x, y); // FIXME
            // case 73:  return iff(x, y); // FIXME

            default:
            msg << "Unknown binary operation: \"" << i << "\"";
          }
          break;
        }

        // N-ary operator, classes 2, 6 and 11 in Gay2005
        case 11: case 12: case 54: case 59: case 60: case 61: case 70: case 71: case 74:
        {
          // Number of elements in the sum
          int n;
          s_ >> n;

          // Collect the arguments
          vector<MX> args(n);
          for (int k=0; k<n; ++k) {
            args[k] = expr();
          }

          // Perform the operation
          switch (i) {
            // case 11: return min(args).scalar(); FIXME // rename?
            // case 12: return max(args).scalar(); FIXME // rename?
            // case 54: return sum(args).scalar(); FIXME // rename?
            // case 59: return count(args).scalar(); FIXME // rename?
            // case 60: return numberof(args).scalar(); FIXME // rename?
            // case 61: return numberofs(args).scalar(); FIXME // rename?
            // case 70: return all(args).scalar(); FIXME // and in AMPL // rename?
            // case 71: return any(args).scalar(); FIXME // or in AMPL // rename?
            // case 74: return alldiff(args).scalar(); FIXME // rename?
            case 54:
            {
              MX r = 0;
              for (vector<MX>::const_iterator it=args.begin();
              it!=args.end(); ++it) r += *it;
              return r;
            }

            default:
            msg << "Unknown n-ary operation: \"" << i << "\"";
          }
          break;
        }

        // Piecewise linear terms, class 4 in Gay2005
        case 64:
        msg << "Piecewise linear terms not supported";
        break;

        // If-then-else expressions, class 5 in Gay2005
        case 35: case 65: case 72:
        msg << "If-then-else expressions not supported";
        break;

        default:
        msg << "Unknown operatio: \"" << i << "\"";
      }
      break;

      default:
      msg << "Unknown instruction: \"" << inst << "\"";
    }

    // Throw error message
    throw CasadiException("Error in NlpBuilder::expr: " + msg.str());
  }

  void NlImporter::F_segment() {
    casadi_error("Imported function description unsupported.");
  }

  void NlImporter::S_segment() {
    casadi_error("Suffix values unsupported");
  }

  void NlImporter::V_segment() {
    // Read header
    int i, j, k;
    s_ >> i >> j >> k;

    // Make sure that v is long enough
    if (i >= v_.size()) {
      v_.resize(i+1);
    }

    // Initialize element to zero
    v_.at(i) = 0;

    // Add the linear terms
    for (int jj=0; jj<j; ++jj) {
      // Linear term
      int pl;
      double cl;
      s_ >> pl >> cl;

      // Add to variable definition (assuming it has already been defined)
      casadi_assert_message(!v_.at(pl).is_empty(), "Circular dependencies not supported");
      v_.at(i) += cl*v_.at(pl);
    }

    // Finally, add the nonlinear term
    v_.at(i) += expr();
  }

  void NlImporter::C_segment() {
    // Get the number
    int i;
    s_ >> i;

    // Parse and save expression
    nlp_.g.at(i) = expr();
  }

  void NlImporter::L_segment() {
    casadi_error("Logical constraint expression unsupported");
  }

  void NlImporter::O_segment() {
    // Get the number
    int i;
    s_ >> i;

    // Should the objective be maximized
    int sigma;
    s_ >> sigma;
    MX sign = sigma!=0 ? -1 : 1;

    // Parse and save expression
    nlp_.f += sign*expr();
  }

  void NlImporter::d_segment() {
    // Read the number of guesses supplied
    int m;
    s_ >> m;

    // Process initial guess for the fual variables
    for (int i=0; i<m; ++i) {
      // Offset and value
      int offset;
      double d;
      s_ >> offset >> d;

      // Save initial guess
      nlp_.lambda_init.at(offset) = d;
    }
  }

  void NlImporter::x_segment() {
    // Read the number of guesses supplied
    int m;
    s_ >> m;

    // Process initial guess
    for (int i=0; i<m; ++i) {
      // Offset and value
      int offset;
      double d;
      s_ >> offset >> d;

      // Save initial guess
      nlp_.x_init.at(offset) = d;
    }
  }

  void NlImporter::r_segment() {
    // For all constraints
    for (int i=0; i<n_con_; ++i) {

      // Read constraint type
      int c_type;
      s_ >> c_type;

      // Temporary
      double c;

      switch (c_type) {
        // Upper and lower bounds
        case 0:
        s_ >> c;
        nlp_.g_lb.at(i) = c;
        s_ >> c;
        nlp_.g_ub.at(i) = c;
        continue;

        // Only upper bounds
        case 1:
        s_ >> c;
        nlp_.g_ub.at(i) = c;
        continue;

        // Only lower bounds
        case 2:
        s_ >> c;
        nlp_.g_lb.at(i) = c;
        continue;

        // No bounds
        case 3:
        continue;

        // Equality constraints
        case 4:
        s_ >> c;
        nlp_.g_lb.at(i) = nlp_.g_ub.at(i) = c;
        continue;

        // Complementary constraints
        case 5:
        {
          // Read the indices
          int ck, ci;
          s_ >> ck >> ci;
          casadi_error("Complementary constraints unsupported");
          continue;
        }

        default:
        casadi_error("Illegal constraint type");
      }
    }
  }

  void NlImporter::b_segment() {
    // For all variable
    for (int i=0; i<n_var_; ++i) {

      // Read constraint type
      int c_type;
      s_ >> c_type;

      // Temporary
      double c;

      switch (c_type) {
        // Upper and lower bounds
        case 0:
        s_ >> c;
        nlp_.x_lb.at(i) = c;
        s_ >> c;
        nlp_.x_ub.at(i) = c;
        continue;

        // Only upper bounds
        case 1:
        s_ >> c;
        nlp_.x_ub.at(i) = c;
        continue;

        // Only lower bounds
        case 2:
        s_ >> c;
        nlp_.x_lb.at(i) = c;
        continue;

        // No bounds
        case 3:
        continue;

        // Equality constraints
        case 4:
        s_ >> c;
        nlp_.x_lb.at(i) = nlp_.x_ub.at(i) = c;
        continue;

        default:
        casadi_error("Illegal variable bound type");
      }
    }
  }

  void NlImporter::k_segment() {
    // Get row offsets
    vector<int> rowind(n_var_+1);

    // Get the number of offsets
    int k;
    s_ >> k;
    casadi_assert(k==n_var_-1);

    // Get the row offsets
    rowind[0]=0;
    for (int i=0; i<k; ++i) {
      s_ >> rowind[i+1];
    }
  }

  void NlImporter::J_segment() {
    // Get constraint number and number of terms
    int i, k;
    s_ >> i >> k;

    // Get terms
    for (int kk=0; kk<k; ++kk) {
      // Get the term
      int j;
      double c;
      s_ >> j >> c;

      // Add to constraints
      nlp_.g.at(i) += c*v_.at(j);
    }
  }

  void NlImporter::G_segment() {
    // Get objective number and number of terms
    int i, k;
    s_ >> i >> k;

    // Get terms
    for (int kk=0; kk<k; ++kk) {
      // Get the term
      int j;
      double c;
      s_ >> j >> c;

      // Add to objective
      nlp_.f += c*v_.at(j);
    }
  }


} // namespace casadi
