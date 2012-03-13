/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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
#include <fstream>

using namespace std;
namespace CasADi{

void SymbolicNLP::parseNL(const std::string& filename, const Dictionary& options){
  // Default options
  bool verbose=false;
  
  // Read user options
  for(Dictionary::const_iterator it=options.begin(); it!=options.end(); ++it){
    if(it->first.compare("verbose")==0){
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
  if(verbose) cout << "Reading file \"" << filename << "\"" << endl;
    
  // Read the header of the NL-file (first 10 lines)
  const int header_sz = 10;
  vector<string> header(header_sz);
  for(int k=0; k<header_sz; ++k){
    getline(nlfile,header[k]);
  }
  
  // Assert that the file is not in binary form
  casadi_assert_message(header.at(0).at(0)=='g',"File could not be read, or file is binary format (currently not supported)");
  
  // Get the number of objectives and constraints
  stringstream ss(header[1]);
  int n_var, n_ineq, n_obj, n_eq, n_lcon;
  ss >> n_var >> n_ineq >> n_obj >> n_eq >> n_lcon;
  
  if(verbose){
    cout << "n_var = " << n_var << ", n_ineq  = " << n_ineq << ", n_obj = " << n_obj << ", n_eq = " << n_eq << ", n_lcon = " << n_lcon << endl;
  }
  
  // Allocate variables
  x = ssym("x",n_var);
  
  // Allocate f and c
  f = SXMatrix::nan(n_obj);
  c = SXMatrix::nan(n_ineq+n_eq);
  
  // Process segments
  while(true){
    
    // Read segment key
    char key;
    nlfile >> key;
  
    // Break if end of file
    if(nlfile.eof()) break;
    
    // Process segments
    switch(key){
      // Imported function description
      case 'F':
	if(verbose) cerr << "Imported function description unsupported: ignored" << endl;
	break;
	
      // Suffix values
      case 'S':
	if(verbose) cerr << "Suffix values unsupported: ignored" << endl;
	break;
	
      // Defined variable definition
      case 'V':
	if(verbose) cerr << "Defined variable definition unsupported: ignored" << endl;
	break;
      
      // Algebraic constraint body
      case 'C':
      {
	// Get the number
	int i;
	nlfile >> i;
	
	// Parse and save expression
	c.at(i) = readExpressionNL(nlfile);
	
	break;
      }
      
      // Logical constraint expression
      case 'L':
	if(verbose) cerr << "Logical constraint expression unsupported: ignored" << endl;
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
	f.at(i) = readExpressionNL(nlfile);
	
	// Negate the expression if we maximize
	if(sigma!=0){
	  f.at(i) = -f.at(i);
	}
	
	break;
      }
      
      // Dual initial guess
      case 'd':
	if(verbose) cerr << "Dual initial guess unsupported: ignored" << endl;
	break;
      
      // Primal initial guess
      case 'x':
	if(verbose) cerr << "Primal initial guess unsupported: ignored" << endl;
	break;
      
      // Bounds on algebraic constraint bodies ("ranges")
      case 'r':
	if(verbose) cerr << "Ranges unsupported: ignored" << endl;
	break;
      
      // Bounds on variable
      case 'b':
	if(verbose) cerr << "Bounds on variable unsupported: ignored" << endl;
	break;
      
      // Jacobian colun counts
      case 'k':
	if(verbose) cerr << "Jacobian colun counts unsupported: ignored" << endl;
	break;
      
      // Jacobian sparsity, linear terms
      case 'J':
	if(verbose) cerr << "Jacobian sparsity, linear terms unsupported: ignored" << endl;
	break;
      
      // Gradient sparsity, linear terms
      case 'G':
	if(verbose) cerr << "Gradient sparsity, linear terms unsupported: ignored" << endl;
	break;
    }
  }
  
  // Close the NL file
  nlfile.close();
}

SX SymbolicNLP::readExpressionNL(std::istream &stream){
  // Read the instruction
  char inst;
  stream >> inst;
  
  // Temporaries
  int i;
  double d;

  // Error message
  stringstream msg;

  // Process instruction
  switch(inst){
    // Symbolic variable
    case 'v':
      // Read the variable number
      stream >> i;
      
      // Return the corresponding expression
      return x.at(i);
    
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
      switch(i){
	
	// Unary operations, class 1 in Gay2005
	case 13:  return floor(readExpressionNL(stream));
	case 14:  return ceil(readExpressionNL(stream));
	case 15:  return abs(readExpressionNL(stream));
	case 16:  return -readExpressionNL(stream);
	case 34:  return !readExpressionNL(stream);
	case 37:  return tanh(readExpressionNL(stream));
	case 38:  return tan(readExpressionNL(stream));
	case 39:  return sqrt(readExpressionNL(stream));
	case 40:  return sinh(readExpressionNL(stream));
	case 41:  return sin(readExpressionNL(stream));
	case 42:  return log10(readExpressionNL(stream));
	case 43:  return log(readExpressionNL(stream));
	case 44:  return exp(readExpressionNL(stream));
	case 45:  return cosh(readExpressionNL(stream));
	case 46:  return cos(readExpressionNL(stream));
	// case 47:  return atanh(readExpressionNL(stream)); FIXME
	case 49:  return atan(readExpressionNL(stream));
	// case 50:  return asinh(readExpressionNL(stream)); FIXME
	case 51:  return asin(readExpressionNL(stream));
	// case 52:  return acosh(readExpressionNL(stream)); FIXME
	case 53:  return acos(readExpressionNL(stream));
	
	// Binary operations, class 2 in Gay2005
	case 0:   return readExpressionNL(stream) + readExpressionNL(stream);
	case 1:   return readExpressionNL(stream) - readExpressionNL(stream);
	case 2:   return readExpressionNL(stream) * readExpressionNL(stream);
	case 3:   return readExpressionNL(stream) / readExpressionNL(stream);
	// case 4:   return rem(readExpressionNL(stream),readExpressionNL(stream)); FIXME
	case 5:   return pow(readExpressionNL(stream),readExpressionNL(stream));
	// case 6:   return readExpressionNL(stream) < readExpressionNL(stream); TODO: Verify this, what is the difference to 'le' == 23 below?
	case 20:  return readExpressionNL(stream) || readExpressionNL(stream);
	case 21:  return readExpressionNL(stream) && readExpressionNL(stream);
	case 22:  return readExpressionNL(stream) < readExpressionNL(stream);
	case 23:  return readExpressionNL(stream) <= readExpressionNL(stream);
	case 24:  return readExpressionNL(stream) == readExpressionNL(stream);
	case 28:  return readExpressionNL(stream) >= readExpressionNL(stream);
	case 29:  return readExpressionNL(stream) > readExpressionNL(stream);
	case 30:  return readExpressionNL(stream) != readExpressionNL(stream);
	// case 48:  return atan2(readExpressionNL(stream),readExpressionNL(stream)); // FIXME
	// case 55:  return intdiv(readExpressionNL(stream),readExpressionNL(stream)); // FIXME
	// case 56:  return precision(readExpressionNL(stream),readExpressionNL(stream)); // FIXME
	// case 57:  return round(readExpressionNL(stream),readExpressionNL(stream)); // FIXME
	// case 58:  return trunc(readExpressionNL(stream),readExpressionNL(stream)); // FIXME
	// case 73:  return iff(readExpressionNL(stream),readExpressionNL(stream)); // FIXME
	
	// N-ary operator, classes 2,6 and 11 in Gay2005
	case 11: case 12: case 54: case 59: case 60: case 61: case 70: case 71: case 74:
	{
	  // Number of elements in the sum
	  int n;
	  stream >> n;
	  
	  // Collect the arguments
	  vector<SX> args(n);
	  for(int k=0; k<n; ++k){
	    args[k] = readExpressionNL(stream);
	  }
	  
	  // Perform the operation
	  switch(i){
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
	      SX r = 0;
	      for(vector<SX>::const_iterator it=args.begin(); it!=args.end(); ++it) r += *it;
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
  
void SymbolicNLP::print(std::ostream &stream) const{
  stream << "NLP:" << endl;
  stream << "x = " << x << endl;
  stream << "#f=" << f.size() << endl;
  stream << "#c=" << c.size() << endl;
}

void SymbolicNLP::repr(std::ostream &stream) const{
  stream << "NLP(#f=" << f.size() << ",#c="<< c.size() << ")";
}



} // namespace CasADi
