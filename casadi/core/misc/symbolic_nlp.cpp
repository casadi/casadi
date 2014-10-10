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


#include "symbolic_nlp.hpp"
#include "../core.hpp"
#include <fstream>

using namespace std;
namespace casadi {

void SymbolicNLP::parseNL(const std::string& filename, const Dictionary& options) {
  // Note: The implementation of this function follows the
  // "Writing .nl Files" paper by David M. Gay (2005)

  // Default options
  bool verbose=false;

  // Read user options
  for (Dictionary::const_iterator it=options.begin(); it!=options.end(); ++it) {
    if (it->first.compare("verbose")==0) {
      verbose = it->second;
    } else {
      stringstream ss;
      ss << "Unknown option \"" << it->first << "\"" << endl;
      throw CasadiException(ss.str());
    }
  }

  // Open the NL file for reading
  ifstream nlfile;
  nlfile.open(filename.c_str());
  if (verbose) cout << "Reading file \"" << filename << "\"" << endl;

  // Read the header of the NL-file (first 10 lines)
  const int header_sz = 10;
  vector<string> header(header_sz);
  for (int k=0; k<header_sz; ++k) {
    getline(nlfile, header[k]);
  }

  // Assert that the file is not in binary form
  casadi_assert_message(header.at(0).at(0)=='g',
                        "File could not be read, or file is binary format "
                        "(currently not supported)");

  // Get the number of objectives and constraints
  stringstream ss(header[1]);
  int n_var, n_con, n_obj, n_eq, n_lcon;
  ss >> n_var >> n_con >> n_obj >> n_eq >> n_lcon;

  if (verbose) {
    cout << "n_var = " << n_var << ", n_con  = " << n_con << ", n_obj = " << n_obj
         << ", n_eq = " << n_eq << ", n_lcon = " << n_lcon << endl;
  }

  // Allocate variables
  x = SX::sym("x", n_var);

  // Allocate f and c
  f = SX::zeros(n_obj);
  g = SX::zeros(n_con);

  // Allocate bounds for x and primal initial guess
  x_lb = DMatrix(x.sparsity(), -numeric_limits<double>::infinity());
  x_ub = DMatrix(x.sparsity(), numeric_limits<double>::infinity());
  x_init = DMatrix(x.sparsity(), 0.0);

  // Allocate bounds for g and dual initial guess
  g_lb = DMatrix(g.sparsity(), -numeric_limits<double>::infinity());
  g_ub = DMatrix(g.sparsity(), numeric_limits<double>::infinity());
  lambda_init = DMatrix(g.sparsity(), 0.0);

  // All variables, including dependent
  vector<SXElement> v = x.data();

  // Process segments
  while (true) {

    // Read segment key
    char key;
    nlfile >> key;

    // Break if end of file
    if (nlfile.eof()) break;

    // Process segments
    switch (key) {
      // Imported function description
      case 'F':
        if (verbose) cerr << "Imported function description unsupported: ignored" << endl;
        break;

      // Suffix values
      case 'S':
        if (verbose) cerr << "Suffix values unsupported: ignored" << endl;
        break;

      // Defined variable definition
      case 'V':
      {
        // Read header
        int i, j, k;
        nlfile >> i >> j >> k;

        // Make sure that v is long enough
        if (i >= v.size()) {
          v.resize(i+1);
        }

        // Initialize element to zero
        v.at(i) = 0;

        // Add the linear terms
        for (int jj=0; jj<j; ++jj) {
          // Linear term
          int pl;
          double cl;
          nlfile >> pl >> cl;

          // Add to variable definition (assuming it has already been defined)
          casadi_assert_message(!v.at(pl).isNan(), "Circular dependencies not supported");
          v[i] += cl*v[pl];
        }

        // Finally, add the nonlinear term
        v[i] += readExpressionNL(nlfile, v);

        break;
      }

      // Algebraic constraint body
      case 'C':
      {
        // Get the number
        int i;
        nlfile >> i;

        // Parse and save expression
        g.at(i) = readExpressionNL(nlfile, v);

        break;
      }

      // Logical constraint expression
      case 'L':
        if (verbose) cerr << "Logical constraint expression unsupported: ignored" << endl;
        break;

      // Objective function
      case 'O':
      {
        // Get the number
        int i;
        nlfile >> i;

        // Should the objective be maximized
        int sigma;
        nlfile >> sigma;

        // Parse and save expression
        f.at(i) = readExpressionNL(nlfile, v);

        // Negate the expression if we maximize
        if (sigma!=0) {
          f.at(i) = -f.at(i);
        }

        break;
      }

      // Dual initial guess
      case 'd':
      {
        // Read the number of guesses supplied
        int m;
        nlfile >> m;

        // Process initial guess for the fual variables
        for (int i=0; i<m; ++i) {
          // Offset and value
          int offset;
          double d;
          nlfile >> offset >> d;

          // Save initial guess
          lambda_init.at(offset) = d;
        }

        break;
      }

      // Primal initial guess
      case 'x':
      {
        // Read the number of guesses supplied
        int m;
        nlfile >> m;

        // Process initial guess
        for (int i=0; i<m; ++i) {
          // Offset and value
          int offset;
          double d;
          nlfile >> offset >> d;

          // Save initial guess
          x_init.at(offset) = d;
        }

        break;
      }

      // Bounds on algebraic constraint bodies ("ranges")
      case 'r':
      {
        // For all constraints
        for (int i=0; i<n_con; ++i) {

          // Read constraint type
          int c_type;
          nlfile >> c_type;

          // Temporary
          double c;

          switch (c_type) {
            // Upper and lower bounds
            case 0:
              nlfile >> c;
              g_lb.at(i) = c;
              nlfile >> c;
              g_ub.at(i) = c;
              continue;

            // Only upper bounds
            case 1:
              nlfile >> c;
              g_ub.at(i) = c;
              continue;

            // Only lower bounds
            case 2:
              nlfile >> c;
              g_lb.at(i) = c;
              continue;

            // No bounds
            case 3:
              continue;

            // Equality constraints
            case 4:
              nlfile >> c;
              g_lb.at(i) = g_ub.at(i) = c;
              continue;

              // Complementary constraints
              case 5:
              {
                // Read the indices
                int ck, ci;
                nlfile >> ck >> ci;
                if (verbose) cerr << "Complementary constraints unsupported: ignored" << endl;
                continue;
              }

              default:
                throw CasadiException("Illegal constraint type");
            }
          }
        break;
      }

      // Bounds on variable
      case 'b':
      {
        // For all variable
        for (int i=0; i<n_var; ++i) {

          // Read constraint type
          int c_type;
          nlfile >> c_type;

          // Temporary
          double c;

          switch (c_type) {
            // Upper and lower bounds
            case 0:
              nlfile >> c;
              x_lb.at(i) = c;
              nlfile >> c;
              x_ub.at(i) = c;
              continue;

            // Only upper bounds
            case 1:
              nlfile >> c;
              x_ub.at(i) = c;
              continue;

           // Only lower bounds
           case 2:
              nlfile >> c;
              x_lb.at(i) = c;
              continue;

           // No bounds
           case 3:
              continue;

           // Equality constraints
           case 4:
              nlfile >> c;
              x_lb.at(i) = x_ub.at(i) = c;
              continue;

           default:
             throw CasadiException("Illegal variable bound type");
          }
        }

        break;
      }

      // Jacobian row counts
      case 'k':
      {
        // Get row offsets
        vector<int> rowind(n_var+1);

        // Get the number of offsets
        int k;
        nlfile >> k;
        casadi_assert(k==n_var-1);

        // Get the row offsets
        rowind[0]=0;
        for (int i=0; i<k; ++i) {
          nlfile >> rowind[i+1];
        }
        break;
      }

      // Linear terms in the constraint function
      case 'J':
      {
        // Get constraint number and number of terms
        int i, k;
        nlfile >> i >> k;

        // Get terms
        for (int kk=0; kk<k; ++kk) {
          // Get the term
          int j;
          double c;
          nlfile >> j >> c;

          // Add to constraints
          g.at(i) += c*v.at(j);
        }
        break;
      }

      // Linear terms in
      case 'G':
      {
        // Get objective number and number of terms
        int i, k;
        nlfile >> i >> k;

        // Get terms
        for (int kk=0; kk<k; ++kk) {
          // Get the term
          int j;
          double c;
          nlfile >> j >> c;

          // Add to objective
          f.at(i) += c*v.at(j);
        }
        break;
      }
    }
  }

  // Close the NL file
  nlfile.close();
}

SXElement SymbolicNLP::readExpressionNL(std::istream &stream, const std::vector<SXElement>& v) {
  // Read the instruction
  char inst;
  stream >> inst;

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
      stream >> i;

      // Return the corresponding expression
      return v.at(i);

    // Numeric expression
    case 'n':

      // Read the floating point number
      stream >> d;

      // Return an expression containing the number
      return d;

    // Operation
    case 'o':

      // Read the operation
      stream >> i;

      // Process
      switch (i) {

        // Unary operations, class 1 in Gay2005
        case 13:  case 14:  case 15:  case 16:  case 34:  case 37:  case 38:  case 39:  case 40:
        case 41:  case 43:  case 42:  case 44:  case 45:  case 46:  case 47:  case 49:  case 50:
        case 51:  case 52:  case 53:
        {
          // Read dependency
          SXElement x = readExpressionNL(stream, v);

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
          SXElement x = readExpressionNL(stream, v);
          SXElement y = readExpressionNL(stream, v);

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
          stream >> n;

          // Collect the arguments
          vector<SXElement> args(n);
          for (int k=0; k<n; ++k) {
            args[k] = readExpressionNL(stream, v);
          }

          // Perform the operation
          switch (i) {
            // case 11: return min(args).toScalar(); FIXME // rename?
            // case 12: return max(args).toScalar(); FIXME // rename?
            // case 54: return sum(args).toScalar(); FIXME // rename?
            // case 59: return count(args).toScalar(); FIXME // rename?
            // case 60: return numberof(args).toScalar(); FIXME // rename?
            // case 61: return numberofs(args).toScalar(); FIXME // rename?
            // case 70: return all(args).toScalar(); FIXME // and in AMPL // rename?
            // case 71: return any(args).toScalar(); FIXME // or in AMPL // rename?
            // case 74: return alldiff(args).toScalar(); FIXME // rename?
            case 54:
            {
              SXElement r = 0;
              for (vector<SXElement>::const_iterator it=args.begin();
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
  throw CasadiException("Error in SymbolicNLP::readExpressionNL: " + msg.str());
}

void SymbolicNLP::print(std::ostream &stream, bool trailing_newline) const {
  stream << "NLP:" << endl;
  stream << "x = " << x << endl;
  stream << "#f=" << f.size() << endl;
  stream << "#g=" << g.size() << endl;
  if (trailing_newline) stream << endl;
}

void SymbolicNLP::repr(std::ostream &stream, bool trailing_newline) const {
  stream << "NLP(#f=" << f.size() << ",#g="<< g.size() << ")";
  if (trailing_newline) stream << endl;
}

} // namespace casadi
