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

#ifndef CASADI_OPTIONS_HPP
#define CASADI_OPTIONS_HPP

namespace CasADi {
  /**
  * \brief Collects global CasADi options
  *
  *
  * Note to developers:  \n
  *  - use sparingly. Global options are - in general - a rather bad idea \n
  *  - this class must never be instantiated. Access its static members directly \n 
  *
  *  \author Joris Gillis 
  *  \date 2012
  */
  class CasadiOptions {
    private:
      /// No instances are allowed
      CasadiOptions();
    public:

      /** \brief Catch CasADi errors when they reach the python layer
      *  If set to true (which will always be the default), any CasADi errors are caught in the SWIG interface
      *  and converted to python exceptions. This allows the user to obtain python stacktraces, use try/except, etc...
      *
      *  If set to false, CasADi errors crash the python interface. This allows a CasADi developer
      *  to obtain a full C++ stacktrace with a debugger such as 'gdb'.
      *
      *  Default: true
      */
      static bool catch_errors_python;
      
      // Setter and getter for catch_errors_python
      static void setCatchErrorsPython(bool flag) { catch_errors_python = flag; }
      static bool getCatchErrorsPython() { return catch_errors_python; }

      /** \brief  Implement the (in)equality operator.
      *  If set to true,  operator==  will check if two casadi objects point to the same node.
      *  If set to false, operator==  will invariably throw an error.
      *
      *  Why this option?
      *  
      *  The mathematician/engineer who wants to use casadi and does not care much about python goodies:
      *  he/she needs to be guided by helpful error messages to avoid doing the incorrect thing.
      *    For example, doing an if-test on a symbolic expression is a common mistake.
      *         'if x==y' with x and y both SX expressions
      *    It's good if an error is thrown in that case.
      *   
      *  ==> Having the option false, is engineer-friendly
      *
      *  However, the ==operator has special meaning in python. It is used by the index method
      *  of lists, used in set maniplations, used in dictionaries.
      *  Having operator== throw an error is very unpythonic
      * 
      *  ==> Having the option true, is python-friendly
      *
      *  Default: false
      */
      static bool equality_operator;
      
      // Setter and getter for equality_operator
      static void setEqualityOperator(bool flag) { equality_operator = flag; }
      static bool getEqualityOperator() { return equality_operator; }
      
  };

}

#endif //CASADI_OPTIONS_HPP
